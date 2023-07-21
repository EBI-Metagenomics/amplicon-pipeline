

import argparse
from collections import defaultdict
from functools import partial
from multiprocessing import Pool
import subprocess

import pandas as pd

from utils import split_dir_into_sample_paths, get_read_count


parser = argparse.ArgumentParser()

parser.add_argument("-d", "--dir", required=True, type=str, help="Input directory with both strands for every sample")
parser.add_argument("-f", "--fwd", required=True, type=str, help="Forward primer sequence")
parser.add_argument("-r", "--rev", required=True, type=str, help="Reverse primer sequence")
parser.add_argument("-o", "--out", required=True, type=str, help="Output directory for cutadapt fastq.gz files")
parser.add_argument("-cpu", "--cpu", type=int, default=1, help="Number of CPUs")

args = parser.parse_args()
_DIR = args.dir
_FWD = args.fwd
_REV = args.rev
_OUT = args.out
_CPU = args.cpu


def cutadapt_one_sample(sample, fwd_primer, rev_primer):
    
# cutadapt -g CCTACGGGNGGCWGCAG -G GACTACHVGGGTATCTAATCC 
# -o cutadapt/SRR13978073_1.cutadapt.fastq -p cutadapt/SRR13978073_2.cutadapt.fastq 
# SRR13978073_1.fastq.gz SRR13978073_2.fastq.gz  

    sample_name = sample.split('/')[-1]

    cmd = [
        'cutadapt',
        '-g',
        fwd_primer,
        '-G',
        rev_primer,
        '-o',
        f'{_OUT}/{sample_name}_1.cutadapt.fastq.gz',
        '-p',
        f'{_OUT}/{sample_name}_2.cutadapt.fastq.gz',
        f'{_DIR}/{sample_name}_1.fastq.gz',
        f'{_DIR}/{sample_name}_2.fastq.gz'
    ]

    print(' '.join(cmd))

    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = proc.communicate()

    print(stdout)

def main():

    sample_list = split_dir_into_sample_paths(_DIR)
    sample_names = [ sample.split('/')[-1] for sample in sample_list ]
    func_partial = partial(cutadapt_one_sample, fwd_primer=_FWD, rev_primer=_REV)

    p = Pool(_CPU) # Prepare pool for parallelisation on per sample basis
    p.map(func_partial, sample_list)

    
if __name__ == "__main__":
    main()