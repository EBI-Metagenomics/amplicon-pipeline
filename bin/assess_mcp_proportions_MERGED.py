
import argparse
from collections import defaultdict
import subprocess

import pandas as pd
import numpy as np

from utils import get_read_count, build_cons_seq

def parse_args():

    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--input", required=True, type=str, help="Path to fastq file to assess mcps")
    parser.add_argument("-s", "--sample", required=True, type=str, help="Sample ID")
    parser.add_argument("-st", "--strand", required=True, choices=['FR', 'F', 'R'], help='F: Forward, R: Reverse')
    parser.add_argument("-o", "--output", required=True, type=str, help="Output path")

    args = parser.parse_args()
    
    _PATH = args.input
    _SAMPLE = args.sample
    _STRAND = args.strand
    _OUTPUT = args.output

    return _PATH, _SAMPLE, _STRAND, _OUTPUT

def fetch_mcp(file, prefix_len, beg=1, rev=False):

    beg = str(beg)
    prefix_len = str(prefix_len)

    if not rev:
        rev = '0'
    else:
        rev = '1'

    cmd = [
        'bash',
        '/hps/software/users/rdf/metagenomics/service-team/users/chrisata/asv_gen/bin/find_most_common_prefixes.sh',
        '-i',
        file,
        '-l',
        prefix_len,
        '-c',
        '10',
        '-b',
        beg,
        '-r',
        rev
    ]

    # print(' '.join(cmd))
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = proc.communicate()
    output = stdout.decode('ascii')
    
    output = output.strip()
    output = output.split('\n')

    mcp_count_dict = defaultdict(int)

    for line in output:
        line = line.strip()
        temp_lst = line.split(' ')
        if len(temp_lst) == 2:
            mcp_count_dict[temp_lst[1]] = int(temp_lst[0])

    return mcp_count_dict


def find_mcp_props_for_sample(_PATH, rev=False):

    res_dict = defaultdict(float)
    lengths_range = range(2, 25, 1)
    
    print(f'Processing {_PATH}')

    subs_len = 4

    for beg in lengths_range:

        l = beg+subs_len
        l = str(l)

        mcp_count_dict = fetch_mcp(_PATH, l, beg, rev)
        mcp_cons_list = []

        mcp_len = len(list(mcp_count_dict.keys())[0])
        for i in range(mcp_len):
            index_base_dict = defaultdict(int)
            for mcp in mcp_count_dict.keys():
                if len(mcp) <= subs_len:
                    continue
                base = mcp[i]
                index_base_dict[base] += mcp_count_dict[mcp]
            mcp_cons_list.append(index_base_dict)

        read_count = get_read_count(_PATH, type='fastq')

        cons_seq, cons_conf = build_cons_seq(mcp_cons_list, read_count)
        
        # mcp_sum = sum(mcp_count_dict.values())

        # fwd_prop = mcp_sum/read_count

        # res_dict_fwd[beg] = fwd_prop
        # res_dict_rev[beg] = rev_prop

        res_dict[beg] = np.mean(cons_conf)

    return res_dict

def concat_out(fwd_out='', rev_out=''):

    total_res_dict = defaultdict(list)
    df_ind = []

    if fwd_out != '':
        [ total_res_dict[key].append(fwd_out[key]) for key in fwd_out.keys() ]
        df_ind.append('F')

    if rev_out != '':
        [ total_res_dict[key].append(rev_out[key]) for key in rev_out.keys() ] 
        df_ind.append('R')

    res_df= pd.DataFrame.from_dict(total_res_dict)

    res_df.index = df_ind

    return res_df


def main():

    _PATH, _SAMPLE, _STRAND, _OUTPUT = parse_args()

    res_df = ''

    match _STRAND:
        case "FR":
            fwd_out = find_mcp_props_for_sample(_PATH)
            rev_out = find_mcp_props_for_sample(_PATH, rev=True)
            res_df = concat_out(fwd_out, rev_out)
            print(res_df)
        case "F":
            fwd_out = find_mcp_props_for_sample(_PATH)
            res_df = concat_out(fwd_out)
            print(res_df)
        case "R":
            rev_out = find_mcp_props_for_sample(_PATH, rev=True)
            res_df = concat_out(rev_out=rev_out)
            print(res_df)
        case _:
            print("Incorrect strand input. Should be F for forward, R for reverse, or FR for both.")
            exit(1)

        
    res_df.to_csv(f'{_OUTPUT}/{_SAMPLE}_mcp_cons.tsv', sep='\t')

    
if __name__ == "__main__":
    main()