
import os
import subprocess

import numpy as np

from constants.regex_ambiguous_bases import _AMBIGUOUS_BASES_DICT

def split_dir_into_sample_paths(_DIR):

    file_list = os.listdir(_DIR)
    file_list = [ file for file in file_list if '.fastq' in file and ('_1' in file or '_2' in file) ]
    sample_set = set()
    [ sample_set.add(f"{_DIR}/{file.split('_')[0]}") for file in file_list ]
    sample_list = sorted(list(sample_set))

    return sample_list

def get_read_count(read_path, type='fastq'):

    cmd = []
    
    if type == 'fastq':
        cmd = [
            'zgrep',
            '-c',
            '^@',
            read_path
        ]

    elif type == 'fasta':
        cmd = [
            'grep',
            '-c',
            '^>',
            read_path
        ]

    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = proc.communicate()

    read_count = int(stdout.strip())

    return read_count

def build_cons_seq(cons_list, read_count, cons_threshold=0.95, do_not_include=[], counter=1):

    cons_seq = ''
    cons_confs = []

    for count_dict in cons_list:
        max_base = '*'
        max_count = 0
        total_count = float(sum(count_dict.values()))
        
        if counter in do_not_include:
            counter += 1
            cons_seq += 'N'
            continue 
        
        for base, count in count_dict.items():
            if base not in ('A', 'T', 'C', 'G'):
                # read_count -= count
                continue

            if count > max_count:
                max_base = base
                max_count = count

        counter += 1

        try:
            prop = max_count/read_count
        except ZeroDivisionError:
            prop = 0.0
        if prop >= cons_threshold:
            cons_seq += max_base
        else:
            cons_seq += 'N'
        cons_confs.append(prop)


    return cons_seq, cons_confs

def primer_regex_query_builder(primer):

    query = ''

    for char in primer:
        if char in ('A', 'C', 'T', 'G'):
            query += char
        else:
            query += str(_AMBIGUOUS_BASES_DICT[char])

    return query
