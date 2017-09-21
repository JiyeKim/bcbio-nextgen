"""low allele fraction variants detection using smCounter

https://github.com/xuchang116/smCounter
"""

import os
import shutil
import sys

from bcbio import utils
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline import datadict as dd
from bcbio.variation import vcfutils
from bcbio.provenance import do

def run(items, background=None):
    print("starting smCounter...")
    return


