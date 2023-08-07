
import argparse
from collections import defaultdict
import subprocess

import pandas as pd
import numpy as np

from utils import get_read_count, build_cons_seq, build_mcp_cons_dict_list

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

def fetch_mcp(fastq, prefix_len, start=1, rev=False):
    """
    Runs a the "find_most_common_prefixes.sh" script to generate the mcps from a fastq file.

    Outputs dictionary containing counts for each generated MCP in the fastq.

    """

    start = str(start)
    prefix_len = str(prefix_len)

    # Check strand requested
    if not rev:
        rev = '0'
    else:
        rev = '1'

    cmd = [
        'bash',
        '/hps/software/users/rdf/metagenomics/service-team/users/chrisata/asv_gen/bin/find_most_common_prefixes.sh',
        '-i',
        fastq,
        '-l',
        prefix_len,
        '-c',
        '10',
        '-b',
        start,
        '-r',
        rev
    ]

    # Run command
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = proc.communicate()
    
    # Capture stdout
    output = stdout.decode('ascii')
    output = output.strip()
    output = output.split('\n')

    mcp_count_dict = defaultdict(int)

    # Process output into mcp count dictionary
    for line in output:
        line = line.strip()
        temp_lst = line.split(' ')
        if len(temp_lst) == 2:
            mcp_count_dict[temp_lst[1]] = int(temp_lst[0])

    return mcp_count_dict


def find_mcp_props_for_sample(_PATH, rev=False):
    """
    Generate mcp proportions in a stepwise and windowed manner for a fastq file.

    For a continuous range of starting indices (2 to 25), generate mcps of window size of 5 bases.
    Calculate the average conservation of the most common base at each index of a window.
    The resulting list of mcp conservations can be considered a conservation curve and used to
    identify inflection points where the conservation suddenly changes.
    
    Output a dictionary where:
        key -> an index starting point e.g. base 10
        val -> the average conservation of the most common base for the mcp window goign from base 10 to 15 (inclusive)
    """

    res_dict = defaultdict(float)
    start_range = range(2, 25, 1) # Range of starting indices
    
    print(f'Processing {_PATH}')

    mcp_len = 5 # length of generated mcps

    for start in start_range:

        end = start+mcp_len-1 # compute the final index for the mcp (inclusive). Indices are of base 1 not 0.
        end = str(end)

        read_count = get_read_count(_PATH, type='fastq') # get read count for fastq file
        mcp_count_dict = fetch_mcp(_PATH, end, start, rev) # get MCP count dict 
        mcp_cons_list = build_mcp_cons_dict_list(mcp_count_dict, mcp_len) # list of base conservation dicts for mcps
        cons_seq, cons_conf = build_cons_seq(mcp_cons_list, read_count) # get list of max base conservations for each index
        
        res_dict[start] = np.mean(cons_conf) # compute the mean

    return res_dict

def concat_out(fwd_out='', rev_out=''):
    """
    Generate Pandas dataframe out of mcp dictionary.

    Output looks like this (when both F and R are requested):
        2	3	4
    F	0.7814975041597337	0.8736772046589019	0.9434276206322796
    R	0.9010981697171381	0.9082861896838601	0.90369384359401

    Columns are the starting indices. Row labels are the strand.
    """

    total_res_dict = defaultdict(list)
    df_ind = []

    # Check if fwd strand was requested
    if fwd_out != '':
        [ total_res_dict[key].append(fwd_out[key]) for key in fwd_out.keys() ]
        df_ind.append('F')

    # Check if rev strand was requested
    if rev_out != '':
        [ total_res_dict[key].append(rev_out[key]) for key in rev_out.keys() ] 
        df_ind.append('R')

    res_df= pd.DataFrame.from_dict(total_res_dict)
    res_df.index = df_ind

    return res_df


def main():

    _PATH, _SAMPLE, _STRAND, _OUTPUT = parse_args()

    res_df = ''

    match _STRAND: # Check for which strands need to be processed
        case "FR": # Both forward and reverse
            fwd_out = find_mcp_props_for_sample(_PATH)
            rev_out = find_mcp_props_for_sample(_PATH, rev=True)
            res_df = concat_out(fwd_out, rev_out)
        case "F": # Only forward
            fwd_out = find_mcp_props_for_sample(_PATH)
            res_df = concat_out(fwd_out)
        case "R": # Only reverse
            rev_out = find_mcp_props_for_sample(_PATH, rev=True)
            res_df = concat_out(rev_out=rev_out)
        case _:
            print("Incorrect strand input. Should be F for forward, R for reverse, or FR for both.")
            exit(1)

    # Save resulting dataframe to a tsv file
    res_df.to_csv(f'{_OUTPUT}/{_SAMPLE}_mcp_cons.tsv', sep='\t')

    
if __name__ == "__main__":
    main()