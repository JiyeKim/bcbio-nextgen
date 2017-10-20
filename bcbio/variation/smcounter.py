"""Barcode-aware variant calling with smCounter.
"""
import os
import sys

import toolz as tz
import pandas as pd

from bcbio import utils
from bcbio.distributed.transaction import file_transaction
from bcbio.heterogeneity import chromhacks
from bcbio.pipeline import config_utils, shared, sample
from bcbio.pipeline import datadict as dd
from bcbio.provenance import do
from bcbio.variation import annotation, bedutils, ploidy, vcfutils
from bcbio.variation.vcfutils import (get_paired_bams, is_paired_analysis,
                                      move_vcf)


def pairwise_alignment(target, real, min_dist=1):
    ''' compare 2 barcodes weather they are same or not.
        return "real" tag if they are same'''
    # print '==================='
    from Bio import pairwise2
    
    if abs(len(target) - len(real)) > min_dist:
        return None
    
    min_score = len(target) - min_dist if len(target) > len(real) else len(real) - min_dist
        
    # print 'Min score:', min_score    
    pairwise2.MAX_ALIGNMENTS = 10
    result = None
    for alignment in pairwise2.align.globalms(target, real, 1, -1, 0, 0):
        # print pairwise2.format_alignment(*alignment)
        if alignment[2] >= min_score:
            result = real
        else:
            continue
    return result

def parse_cigars_get_barcode_pos(cigar, flag):
    import re
    cigar_pattern = '([0-9]+[MIDNSHP])'
    cigar_search = re.compile(cigar_pattern)
    cigar_strings = cigar_search.findall(cigar)
    atomic_cigar_pattern = '([0-9]+)([MIDNSHP])'
    atomic_cigar_search = re.compile(atomic_cigar_pattern)

    if flag == 0: # +strand
        cigar_strings.reverse()
        
    barcode_pos = 0
    for cigar_string in cigar_strings:
        cigar_op_list = atomic_cigar_search.findall(cigar_string)[0]
        if cigar_op_list[1] == 'M':
            break
        else:
            barcode_pos += int(cigar_op_list[0])
            
    return (cigar_strings, barcode_pos)


def get_12n_barcode(row, n=12):
    # print row.name
    if row.FLAG not in [0, 16]:
        return None
    
    strings, barcode_pos = parse_cigars_get_barcode_pos(row.CIGAR, row.FLAG)
    
    # barcode is located upstream when read is - strand
    if barcode_pos == 0:
        return None
    barcode = row.SEQ[:barcode_pos] if row.FLAG == 16 else row.SEQ[-barcode_pos:]
    
    # cut barcode sequence 12 base
    if len(barcode) > 12:
        if row.FLAG == 0: # + strand
            barcode = barcode[:n]
        elif row.FLAG == 16: # - strand
            barcode = barcode[-n:]
            
    # check barcode length (default: > 7bp)
    minimum_barcode_length = n - 5
    if len(barcode) < minimum_barcode_length:
        return None
    return barcode


def compare_and_merge_barcode(target, barcodes_df, min_dist):
    ''' return matched barcodes in realgroup with given target '''
    if target.is_real is True:
        return target.name
    elif target.modified is not None:
        return target.modified
    
    modified = None
    for realb in barcodes_df[barcodes_df.is_real == True].itertuples():
        modified = pairwise_alignment(target.name, realb.Index, min_dist)
        if modified:
            ''' check the number of reads is large enough to merge '''
            if target.reads_no * 6 <= realb.reads_no:
                total_reads = target.reads_no + realb.reads_no
                barcodes_df.set_value(realb.Index, 'reads_no', total_reads)
                barcodes_df.set_value(target.name, 'reads_no', 0)
                break
            else:
                modified = None
    return modified


def get_new_barcode(row):
    original_barcode = row.BARCODE if row.BARCODE else 'undefined'
    clustered_barcode = row.FINAL_BARCODE if row.FINAL_BARCODE else original_barcode
    
    new_id = '{ori_id}:{chrno}-{strand}-{position}-{cluster_bar}:{ori_bar}'.format(
        ori_id=row.name,
        chrno=row.RNAME,
        strand='1' if row.FLAG == 16 else '0',
        position=row.POS,
        cluster_bar=clustered_barcode,
        ori_bar=original_barcode,
    )
    print(new_id)
    return new_id


def get_final_barcode(row, final_barcodes):
    if row.BARCODE is None:
        return None
    return final_barcodes.loc[row.BARCODE].modified


def _run_cluster_barcode(bamfile):
    prefix = os.path.splitext(bamfile)[0]
    samfile = "{}.sam".format(prefix)

    # convert bam to sam
    do.run("samtools view -G 4 {bamfile} > {samfile}".format(**locals()))
    
    # barcode extraction
    df = pd.read_table(
        samfile, index_col=0, header=None,
        usecols=list(range(0, 11)),
        names = ['QNAME', 'FLAG', 'RNAME', 'POS',
                 'MAPQ', 'CIGAR', 'RNEXT', 'PNEXT',
                 'TLEN', 'SEQ', 'QUAL'],
    )
    df['BARCODE'] = df.apply(get_12n_barcode, axis=1)
    df_dropna = df.dropna(subset=['BARCODE'])

    print 'total reads:', df.shape
    print 'total reads after dropna:', df_dropna.shape

    # grouping
    read_groups = df_dropna.sort_values(
        ['FLAG'], ascending=False).groupby(by=['BARCODE'])
    count = read_groups.agg('count')
    count = count.query('BARCODE != 0')
    count.to_csv("{prefix}_read_groups_sample_by_barcode.csv".format(**locals()))
    print("all groups:", count.index.size)
    print("1 read group:", count.query('FLAG == 1').index.size)
    print("over 1 read group:", count.query('FLAG > 1').index.size)
    print count.FLAG.describe()

    # Barcode clustering
    # step 1
    print(' --- step1 --- ')
    ''' sorting by number of reads in groups '''
    barcodes = count.sort_values(by='FLAG', ascending=False)[['FLAG']]
    barcodes.columns = ['reads_no']
    real_barcode_cutoff = int(barcodes.reads_no.max() * 0.05)
    barcodes['is_real'] = barcodes.reads_no >= real_barcode_cutoff
    barcodes.head()
    print 'real barcodes (reads >', real_barcode_cutoff, "):", \
          barcodes.is_real.sum(), '/', barcodes.index.size

    # Step 2
    print(' --- step2 --- ')
    barcodes['modified'] = None
    barcodes['modified'] = barcodes.apply(
        compare_and_merge_barcode, axis=1, barcodes_df=barcodes, min_dist=1)
    barcodes[:10]
    print 'total real barcodes:', barcodes.modified.count(), \
          '/', barcodes.index.size

    # Step 3
    print(' --- step3 --- ')
    barcodes_step3 = barcodes.copy()
    barcodes_step3.loc[
        (barcodes_step3.is_real == False) & \
        (barcodes_step3.reads_no >= 2), 'is_real'] = True
    print barcodes_step3.head()

    # Step 4
    print(' --- step4 --- ')
    barcodes_step4 = barcodes_step3.copy()
    barcodes_step4['modified'] = barcodes_step4.apply(
        compare_and_merge_barcode, axis=1,
        barcodes_df=barcodes_step4, min_dist=2
    )
    barcodes_step4[:10]
    print 'total real barcodes:', barcodes_step4.modified.count(), \
          '/', barcodes_step4.index.size

    # Summary
    a = barcodes_step4.query('is_real == False and reads_no > 0').index.size
    print('undefined barcodes: ', a)
    b = barcodes_step4.query('is_real == True').index.size
    print('real barcodes:', b)
    c = barcodes_step4.query('is_real == False and reads_no == 0').index.size
    print('merged barcodes:', c)
    print('total:', a + b + c)

    final_barcodes = barcodes_step4.copy()
    df['FINAL_BARCODE'] = df.apply(
        get_final_barcode, axis=1, final_barcodes=final_barcodes)
    df['FINAL_QNAME'] = df.apply(get_new_barcode, axis=1)
    
    print 'writing new samfile...'
    tempsamfile = "{}_temp.sam".format(prefix)
    finalsamfile = "{}_final.sam".format(prefix)
    finalbamfile = "{}_final.bam".format(prefix)
    finalbaifile = "{}_final.bai".format(prefix)
    with open(tempsamfile, 'w') as f:
        for line in open(samfile):
            items = line.split('\t')
            original_id = items[0]
            if not isinstance(df.loc[original_id], pd.Series):
                continue

            new_id = df.loc[original_id].FINAL_QNAME
            line = line.replace(original_id, new_id)
            f.write(line)

    do.run(
        'samtools view -H {bamfile} > {finalsamfile} && \
        cat {tempsamfile} >> {finalsamfile}'.format(**locals())
    )
    do.run('samtools view -Sb {finalsamfile} > {finalbamfile}'.format(**locals()))
    do.run('samtools index {finalbamfile} {finalbaifile}'.format(**locals()))
    print('Done--->', finalbamfile)
    return finalbamfile


def run_smcounter(items):
    # print('in run_smcounter()...')
    # print('items--->>', items)
    align_bams = [dd.get_align_bam(x) for x in [items]]

    # barcode clustering
    align_bams = [_run_cluster_barcode(bamfile) for bamfile in align_bams]

    # print 'align_bams--->', align_bams
    out_file = os.path.join(
        dd.get_work_dir(items), "smcounter", "variantcall.vcf.gz")
    ref_file = dd.get_ref_file(items)
    assoc_files = tz.get_in(("genome_resources", "variation"), items, {})
    if not assoc_files: assoc_files = {}
    region = None

    if is_paired_analysis(align_bams, items):
        paired = get_paired_bams(align_bams, items)
        print("is_paired_analysis--->paired")
        print("Not supported...")
        '''
        if not paired.normal_bam:
            call_file = _run_smcounter_caller(align_bams, items, ref_file,
                                              assoc_files, region, out_file, somatic=paired)
        else:
            call_file = _run_smcounter_paired([paired.tumor_bam, paired.normal_bam],
                                              [paired.tumor_data, paired.normal_data],
                                              ref_file, assoc_files, region, out_file)
        '''

    else:
        print("is_paired_analysis--->single")
        vcfutils.check_paired_problems(items)
        print('outfile----->', out_file)
        call_file = _run_smcounter_caller(align_bams, items, ref_file,
                                          assoc_files, region, out_file)

    return call_file


def _run_smcounter_caller(align_bams, items, ref_file, assoc_files,
                          region=None, out_file=None):
    print("New! function!")
    config = items["config"]
    smcounter = config['resources']['smcounter']['dir']
    input_bams = ' '.join(align_bams)
    bedtarget = config['algorithm']['variant_regions']
    ref_file = ref_file
    bedtoolspath = '/incogwas/igsite/storage/bcbio/anaconda/bin/'
    runpath = os.path.join(items['dirs']['work'], 'variantcall', items["description"])

    # make directory for smCounter results
    if not os.path.exists(runpath):
        os.makedirs(runpath)

    cmd = ("python {smcounter} --outPrefix smcounter --bamFile {input_bams} "
           "--bedTarget {bedtarget} --minBQ 20 --minMQ 30 --hpLen 10 "
           "--rpb 1.5 --mismatchThr 6.0 --mtDrop 0 --primerDist 2 --mtDepth 3 "
           "--refGenome {ref_file} --bedtoolsPath {bedtoolspath} "
           "--runPath {runpath} --logFile {runpath}/smCounterLog")
    cmdstring = cmd.format(**locals())
    # print('cmd----------->', cmdstring)

    do.run(cmdstring, "Genotyping with smCounter", {})
    return out_file
    
def _run_smcounter_paired(align_bams, items, ref_file, assoc_files,
                          region=None, out_file=None):
    """Detect SNPs and indels with smCounter for paired tumor/normal samples.

    Sources of options for smCounter:
    mailing list: https://groups.google.com/d/msg/smcounter/dTWBtLyM4Vs/HAK_ZhJHguMJ
    mailing list: https://groups.google.com/forum/#!msg/smcounter/LLH7ZfZlVNs/63FdD31rrfEJ
    speedseq: https://github.com/cc2qe/speedseq/blob/e6729aa2589eca4e3a946f398c1a2bdc15a7300d/bin/speedseq#L916
    sga/smcounter: https://github.com/jts/sga-extra/blob/7e28caf71e8107b697f9be7162050e4fa259694b/
                   sga_generate_varcall_makefile.pl#L299
    """
    config = items[0]["config"]
    if out_file is None:
        out_file = "%s-paired-variants.vcf.gz" % os.path.splitext(align_bams[0])[0]
    if not utils.file_exists(out_file):
        raw_file = "%s-raw%s" % utils.splitext_plus(out_file)
        if not utils.file_exists(raw_file):
            with file_transaction(items[0], raw_file) as tx_out_file:
                paired = get_paired_bams(align_bams, items)
                assert paired.normal_bam, "Require normal BAM for smCounter paired calling and filtering"

                smcounter = config_utils.get_program("smcounter", config)
                opts, no_target_regions = _smcounter_options_from_config(items, config, out_file, region)
                if no_target_regions:
                    vcfutils.write_empty_vcf(tx_out_file, config,
                                            samples=[x for x in [paired.tumor_name, paired.normal_name] if x])
                else:
                    opts = " ".join(opts)
                    opts += " --min-repeat-entropy 1"
                    opts += " --no-partial-observations"
                    opts = _add_somatic_opts(opts, paired)
                    compress_cmd = "| bgzip -c" if out_file.endswith("gz") else ""
                    # For multi-sample outputs, ensure consistent order
                    samples = ("-s " + ",".join([dd.get_sample_name(d) for d in items])) if len(items) > 1 else ""
                    fix_ambig = vcfutils.fix_ambiguous_cl()
                    bcbio_py = sys.executable
                    py_cl = os.path.join(os.path.dirname(sys.executable), "py")
                    cl = ("{smcounter} -f {ref_file} {opts} "
                          "{paired.tumor_bam} {paired.normal_bam} "
                          """| bcftools filter -i 'ALT="<*>" || QUAL > 5' """
                          """| {bcbio_py} -c 'from bcbio.variation import smcounter; """
                          """smcounter.call_somatic("{paired.tumor_name}", "{paired.normal_name}")' """
                          "| {fix_ambig} | bcftools view {samples} -a - | "
                          "{py_cl} -x 'bcbio.variation.smcounter.remove_missingalt(x)' | "
                          "vcfallelicprimitives -t DECOMPOSED --keep-geno | vcffixup - | vcfstreamsort | "
                          "vt normalize -n -r {ref_file} -q - | vcfuniqalleles | vt uniq - 2> /dev/null "
                          "{compress_cmd} > {tx_out_file}")
                    do.run(cl.format(**locals()), "Genotyping paired variants with smCounter", {})
            ann_file = annotation.annotate_nongatk_vcf(raw_file, align_bams,
                                                       assoc_files.get("dbsnp"),
                                                       ref_file, items[0], out_file)
            if ann_file != out_file:
                utils.symlink_plus(ann_file, out_file)
    return out_file

# ## Filtering

def _check_lods(parts, tumor_thresh, normal_thresh, indexes):
    """Ensure likelihoods for tumor and normal pass thresholds.

    Skipped if no smCounter GL annotations available.
    """
    try:
        gl_index = parts[8].split(":").index("GL")
    except ValueError:
        return True
    try:
        tumor_gls = [float(x) for x in parts[indexes["tumor"]].strip().split(":")[gl_index].split(",") if x != "."]
        if tumor_gls:
            tumor_lod = max(tumor_gls[i] - tumor_gls[0] for i in range(1, len(tumor_gls)))
        else:
            tumor_lod = -1.0
    # No GL information, no tumor call (so fail it)
    except IndexError:
        tumor_lod = -1.0
    try:
        normal_gls = [float(x) for x in parts[indexes["normal"]].strip().split(":")[gl_index].split(",") if x != "."]
        if normal_gls:
            normal_lod = min(normal_gls[0] - normal_gls[i] for i in range(1, len(normal_gls)))
        else:
            normal_lod = normal_thresh
    # No GL inofmration, no normal call (so pass it)
    except IndexError:
        normal_lod = normal_thresh
    return normal_lod >= normal_thresh and tumor_lod >= tumor_thresh

def _check_freqs(parts, indexes):
    """Ensure frequency of tumor to normal passes a reasonable threshold.

    Avoids calling low frequency tumors also present at low frequency in normals,
    which indicates a contamination or persistent error.
    """
    thresh_ratio = 2.7
    try:  # smCounter
        ao_index = parts[8].split(":").index("AO")
        ro_index = parts[8].split(":").index("RO")
    except ValueError:
        ao_index, ro_index = None, None
    try:  # VarDict
        af_index = parts[8].split(":").index("AF")
    except ValueError:
        af_index = None
    if af_index is None and ao_index is None:
        # okay to skip if a gVCF record
        if parts[4].find("<*>") == -1:
            raise NotImplementedError("Unexpected format annotations: %s" % parts[8])
    def _calc_freq(item):
        try:
            if ao_index is not None and ro_index is not None:
                ao = sum([int(x) for x in item.split(":")[ao_index].split(",")])
                ro = int(item.split(":")[ro_index])
                freq = ao / float(ao + ro)
            elif af_index is not None:
                freq = float(item.split(":")[af_index])
            else:
                freq = 0.0
        except (IndexError, ValueError, ZeroDivisionError):
            freq = 0.0
        return freq
    tumor_freq, normal_freq = _calc_freq(parts[indexes["tumor"]]), _calc_freq(parts[indexes["normal"]])
    return normal_freq <= 0.001 or normal_freq <= tumor_freq / thresh_ratio

def remove_missingalt(line):
    """Remove lines that are missing an alternative allele.

    During cleanup of extra alleles, bcftools has an issue in complicated cases
    with duplicate alleles and will end up stripping all alternative alleles.
    This removes those lines to avoid issues downstream.
    """
    if not line.startswith("#"):
        parts = line.split("\t")
        if parts[4] == ".":
            return None
    return line

def call_somatic(tumor_name, normal_name):
    """Call SOMATIC variants from tumor/normal calls, adding REJECT filters and SOMATIC flag.

    Works from stdin and writes to stdout, finding positions of tumor and normal samples.

    Uses MuTect like somatic filter based on implementation in speedseq:
    https://github.com/cc2qe/speedseq/blob/e6729aa2589eca4e3a946f398c1a2bdc15a7300d/bin/speedseq#L62

    Extracts the genotype likelihoods (GLs) from smCounter, which are like phred scores
    except not multiplied by 10.0 (https://en.wikipedia.org/wiki/Phred_quality_score).
    For tumors, we retrieve the best likelihood to not be reference (the first GL) and
    for normal, the best likelhood to be reference.

    After calculating the likelihoods, we compare these to thresholds to pass variants
    at tuned sensitivity/precision. Tuning done on DREAM synthetic 3 dataset evaluations.

    We also check that the frequency of the tumor exceeds the frequency of the normal by
    a threshold to avoid calls that are low frequency in both tumor and normal. This supports
    both smCounter and VarDict output frequencies.
    """
    # Thresholds are like phred scores, so 3.5 = phred35
    tumor_thresh, normal_thresh = 3.5, 3.5
    new_headers = ['##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description="Somatic event">\n',
                   ('##FILTER=<ID=REJECT,Description="Not somatic due to normal call frequency '
                    'or phred likelihoods: tumor: %s, normal %s.">\n')
                   % (int(tumor_thresh * 10), int(normal_thresh * 10))]
    def _output_filter_line(line, indexes):
        parts = line.split("\t")
        if _check_lods(parts, tumor_thresh, normal_thresh, indexes) and _check_freqs(parts, indexes):
            parts[7] = parts[7] + ";SOMATIC"
        else:
            if parts[6] in set([".", "PASS"]):
                parts[6] = "REJECT"
            else:
                parts[6] += ";REJECT"
        line = "\t".join(parts)
        sys.stdout.write(line)
    def _write_header(header):
        for hline in header[:-1] + new_headers + [header[-1]]:
            sys.stdout.write(hline)
    header = []
    indexes = None
    for line in sys.stdin:
        if not indexes:
            if line.startswith("#"):
                header.append(line)
            else:
                parts = header[-1].rstrip().split("\t")
                indexes = {"tumor": parts.index(tumor_name), "normal": parts.index(normal_name)}
                _write_header(header)
                _output_filter_line(line, indexes)
        else:
            _output_filter_line(line, indexes)
    # no calls, only output the header
    if not indexes:
        _write_header(header)

def _clean_smcounter_output(line):
    """Clean smCounter output to make post-processing with GATK happy.

    XXX Not applied on recent versions which fix issues to be more compatible
    with bgzip output, but retained in case of need.

    - Remove lines from smCounter outputs where REF/ALT are identical:
      2       22816178        .       G       G       0.0339196
      or there are multiple duplicate alleles:
      4       60594753        .       TGAAA   T,T
    - Remove Type=Int specifications which are not valid VCF and GATK chokes
      on.
    """
    if line.startswith("#"):
        line = line.replace("Type=Int,D", "Type=Integer,D")
        return line
    else:
        parts = line.split("\t")
        alleles = [x.strip() for x in parts[4].split(",")] + [parts[3].strip()]
        if len(alleles) == len(set(alleles)):
            return line
    return None

def clean_vcf_output(orig_file, clean_fn, config, name="clean"):
    """Provide framework to clean a file in-place, with the specified clean
    function.
    """
    base, ext = utils.splitext_plus(orig_file)
    out_file = "{0}-{1}{2}".format(base, name, ext)
    if not utils.file_exists(out_file):
        with open(orig_file) as in_handle:
            with file_transaction(config, out_file) as tx_out_file:
                with open(tx_out_file, "w") as out_handle:
                    for line in in_handle:
                        update_line = clean_fn(line)
                        if update_line:
                            out_handle.write(update_line)
        move_vcf(orig_file, "{0}.orig".format(orig_file))
        move_vcf(out_file, orig_file)
        with open(out_file, "w") as out_handle:
            out_handle.write("Moved to {0}".format(orig_file))
