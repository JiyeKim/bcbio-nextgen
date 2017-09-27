"""Barcode-aware variant calling with smCounter.

"""

import os
import sys

import toolz as tz

from bcbio import utils
from bcbio.distributed.transaction import file_transaction
from bcbio.heterogeneity import chromhacks
from bcbio.pipeline import config_utils, shared, sample
from bcbio.pipeline import datadict as dd
from bcbio.provenance import do
from bcbio.variation import annotation, bedutils, ploidy, vcfutils
from bcbio.variation.vcfutils import (get_paired_bams, is_paired_analysis,
                                      move_vcf)

def region_to_smcounter(region):
    if isinstance(region, (list, tuple)):
        chrom, start, end = region
        return "%s:%s..%s" % (chrom, start, end)
    else:
        return region

def _smcounter_options_from_config(items, config, out_file, region=None):
    """Prepare standard options from configuration input.

    Input BED target files are merged to avoid overlapping regions which
    cause smCounter to call multiple times.

    Checks for empty sets of target regions after filtering for high depth,
    in which case we should skip the smCounter run.
    """
    opts = ["--genotype-qualities", "--strict-vcf"]
    cur_ploidy = ploidy.get_ploidy(items, region)
    base_ploidy = ploidy.get_ploidy(items)
    opts += ["--ploidy", str(cur_ploidy)]
    # Adjust min fraction when trying to call more sensitively in certain
    # regions. This is primarily meant for pooled mitochondrial calling.
    if (isinstance(region, (list, tuple)) and chromhacks.is_mitochondrial(region[0])
          and cur_ploidy >= base_ploidy and "--min-alternate-fraction" not in opts and "-F" not in opts):
        opts += ["--min-alternate-fraction", "0.01"]
    variant_regions = bedutils.merge_overlaps(bedutils.population_variant_regions(items), items[0])
    # Produce gVCF output
    if any("gvcf" in dd.get_tools_on(d) for d in items):
        opts += ["--gvcf", "--gvcf-chunk", "50000"]
    no_target_regions = False
    target = shared.subset_variant_regions(variant_regions, region, out_file, items)
    if target:
        if isinstance(target, basestring) and os.path.isfile(target):
            if any(tz.get_in(["config", "algorithm", "coverage_interval"], x, "").lower() == "genome"
                   for x in items):
                target = shared.remove_highdepth_regions(target, items)
                if os.path.getsize(target) == 0:
                    no_target_regions = True
            opts += ["--targets", target]
        else:
            opts += ["--region", region_to_smcounter(target)]
    resources = config_utils.get_resources("smcounter", config)
    print('resources----->', resources)
    if resources.get("options"):
        opts += resources["options"]
    return opts, no_target_regions

def _add_somatic_opts(opts, paired):
    """Add somatic options to current set. See _run_smcounter_paired for references.
    """
    if "--min-alternate-fraction" not in opts and "-F" not in opts:
        # add minimum reportable allele frequency
        # smCounter defaults to 20%, but use 10% by default for the
        # tumor case
        min_af = float(utils.get_in(paired.tumor_config, ("algorithm",
                                                          "min_allele_fraction"), 10)) / 100.0
        opts += " --min-alternate-fraction %s" % min_af
    # Recommended settings for cancer calling
    opts += (" --pooled-discrete --pooled-continuous "
             "--report-genotype-likelihood-max --allele-balance-priors-off")
    return opts

def run_smcounter(items):
    print('in run_smcounter()...')
    print('items--->>', items)
    align_bams = [dd.get_align_bam(x) for x in [items]]
    out_file = os.path.join(
        dd.get_work_dir(items), "smcounter", "variantcall.vcf.gz")
    ref_file = dd.get_ref_file(items)
    assoc_files = tz.get_in(("genome_resources", "variation"), items, {})
    if not assoc_files: assoc_files = {}
    region = None

    if is_paired_analysis(align_bams, items):
        paired = get_paired_bams(align_bams, items)
        print("is_paired_analysis--->paired")
        if not paired.normal_bam:
            call_file = _run_smcounter_caller(align_bams, items, ref_file,
                                              assoc_files, region, out_file, somatic=paired)
        else:
            call_file = _run_smcounter_paired([paired.tumor_bam, paired.normal_bam],
                                              [paired.tumor_data, paired.normal_data],
                                              ref_file, assoc_files, region, out_file)
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
    runpath = items['dirs']['work']

    cmd = ("python {smcounter} --outPrefix smcounter --bamFile {input_bams} "
           "--bedTarget {bedtarget} --minBQ 20 --minMQ 30 --hpLen 10 "
           "--rpb 1.5 --mismatchThr 6.0 --mtDrop 0 --primerDist 2 --mtDepth 3 "
           "--refGenome {ref_file} --bedtoolsPath {bedtoolspath} "
       	   "--runPath {runpath} --logFile {runpath}/log")
    cmdstring = cmd.format(**locals())
    print('cmd----------->', cmdstring)

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
