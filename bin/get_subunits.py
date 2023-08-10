#!/usr/bin/env python3

import argparse
import os
from Bio import SeqIO


SSU = "SSU_rRNA"
LSU = "LSU_rRNA"
Seq5S = "mtPerm-5S"
Seq5_8S = "5_8S_rRNA"

SSU_rRNA_archaea = "SSU_rRNA_archaea"
SSU_rRNA_bacteria = "SSU_rRNA_bacteria"
SSU_rRNA_eukarya = "SSU_rRNA_eukarya"
SSU_rRNA_microsporidia = "SSU_rRNA_microsporidia"

LSU_rRNA_archaea = "LSU_rRNA_archaea"
LSU_rRNA_bacteria = "LSU_rRNA_bacteria"
LSU_rRNA_eukarya = "LSU_rRNA_eukarya"


def set_model_names(prefix, name):
    pattern_dict = {}
    pattern_dict[SSU] = os.path.join(directory, f'{name}_SSU.fasta')
    pattern_dict[SSU_rRNA_archaea] = os.path.join(directory, f'{prefix}{name}_{SSU_rRNA_archaea}.RF01959.fa')
    pattern_dict[SSU_rRNA_bacteria] = os.path.join(directory, f'{prefix}{name}_{SSU_rRNA_bacteria}.RF00177.fa')
    pattern_dict[SSU_rRNA_eukarya] = os.path.join(directory, f'{prefix}{name}_{SSU_rRNA_eukarya}.RF01960.fa')
    pattern_dict[SSU_rRNA_microsporidia] = os.path.join(directory, f'{prefix}{name}_{SSU_rRNA_microsporidia}.RF02542.fa')
    pattern_dict[LSU] = os.path.join(directory, f'{name}_LSU.fasta')
    pattern_dict[LSU_rRNA_archaea] = os.path.join(directory, f'{prefix}{name}_{LSU_rRNA_archaea}.RF02540.fa')
    pattern_dict[LSU_rRNA_bacteria] = os.path.join(directory, f'{prefix}{name}_{LSU_rRNA_bacteria}.RF02541.fa')
    pattern_dict[LSU_rRNA_eukarya] = os.path.join(directory, f'{prefix}{name}_{LSU_rRNA_eukarya}.RF02543.fa')
    pattern_dict[Seq5S] = os.path.join(directory, f'{name}_5S.fa')
    pattern_dict[Seq5_8S] = os.path.join(directory, f'{name}_5_8S.fa')
    return pattern_dict


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Extract lsu, ssu and 5s and other models"
    )
    parser.add_argument(
        "-i", "--input", dest="input", help="Input fasta file", required=True
    )
    parser.add_argument(
        "-p", "--prefix", dest="prefix", help="prefix for models", required=False
    )
    parser.add_argument("-n", "--name", dest="name", help="Accession", required=True)

    args = parser.parse_args()
    prefix = args.prefix if args.prefix else ""
    name = args.name if args.name else "accession"

    directory = "sequence-categorisation"
    if not os.path.exists(directory):
        os.makedirs(directory)
    directory_ncRNA = os.path.join("sequence-categorisation", "ncRNA")
    if not os.path.exists(directory_ncRNA):
        os.makedirs(directory_ncRNA)

    print('Start fasta mode')
    pattern_dict = set_model_names(prefix, name)
    coding_rna = [SSU_rRNA_archaea, SSU_rRNA_bacteria, SSU_rRNA_eukarya, SSU_rRNA_microsporidia,
                  LSU_rRNA_archaea, LSU_rRNA_bacteria, LSU_rRNA_eukarya, Seq5S, Seq5_8S]
    open_files = {}
    for record in SeqIO.parse(args.input, "fasta"):
        model = '-'.join(record.id.split('/')[0].split('-')[1:])
        if model in coding_rna:
            filename = pattern_dict[model]
        else:
            filename = os.path.join(directory_ncRNA, f'{prefix}{name}_{model}.fasta')
        if model not in open_files:
            file_out = open(filename, 'w')
            open_files[model] = file_out
        SeqIO.write(record, open_files[model], "fasta")

        if model in (SSU_rRNA_archaea, SSU_rRNA_bacteria, SSU_rRNA_eukarya, SSU_rRNA_microsporidia):
            if SSU not in open_files:
                file_out = open(pattern_dict[SSU], 'w')
                open_files[SSU] = file_out
            SeqIO.write(record, open_files[SSU], "fasta")
        if model in (LSU_rRNA_archaea, LSU_rRNA_bacteria, LSU_rRNA_eukarya):
            if LSU not in open_files:
                file_out = open(pattern_dict[LSU], 'w')
                open_files[LSU] = file_out
            SeqIO.write(record, open_files[LSU], "fasta")

    for item in open_files:
        open_files[item].close()

    if len(os.listdir(directory_ncRNA)) == 0:
        os.rmdir(directory_ncRNA)