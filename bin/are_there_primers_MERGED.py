
import argparse
from collections import defaultdict

import numpy as np

from assess_mcp_proportions_MERGED import fetch_mcp
from utils import get_read_count, build_cons_seq, build_mcp_cons_dict_list

def parse_args():

    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--input", required=True, type=str, help="Path to fastq file to check for primers")
    parser.add_argument("-s", "--sample", required=True, type=str, help="Sample ID")
    parser.add_argument("-o", "--output", required=True, type=str, help="Output path")
    args = parser.parse_args()
    
    _PATH = args.input
    _SAMPLE = args.sample
    _OUTPUT = args.output

    return _PATH, _SAMPLE, _OUTPUT

def are_there_primers_in_this_sample(_PATH, rev=False):
    """
    Predict the presence of primers based on windows of base conservation.

    Takes a fastq file as input. Extracts proportion of most common base for the first 100 bases.
    Computes the a threshold (Q3 - 0.15) based on this proportion and counts the number of bases below 
    it in windows of 10 bases.
    If at least one of the first two windows contains at most one such a base, then the presence of a primer is flagged as true.
    A primer is also flagged as true if the combined count of bases below Q3 is at most 4.

    The output of this function is a boolean flag:
        True if a primer was identified
        False if a primer was not identified
    """

    read_count = get_read_count(_PATH, 'fastq') # Get read count for fastq file
    mcp_len = 100 # Script will look at first 100 base mcps (for rev=True, it will look at first 100 from 3' to 5')

    mcp_count_dict = fetch_mcp(_PATH, mcp_len, rev=rev) # mcp dict where key is the mcp and value is the count
    mcp_cons_list = build_mcp_cons_dict_list(mcp_count_dict, mcp_len) # list of base conservation dicts for mcps
    cons_seq, cons_confs = build_cons_seq(mcp_cons_list, read_count) # get list of max base conservations for each index

  
    window_size = 10
    # Counter that will reset to 0 every 10 bases
    window_count = 0
    # Will append the window count to this list every 10 bases
    window_count_list = []
    # Compute Q3-based threshold
    max_cons = np.quantile(cons_confs, 0.75)
    threshold = max_cons - 0.15 

    if max_cons < 0.75:
        threshold = 0.75
    # Immediately return false (no primer) if the max conservation is less than 0.6
    if max_cons < 0.6:
        return False

    # Loop through every base
    for i, val in enumerate(cons_confs):
        if i%window_size == 0 and i !=0: # After looping through a window..
            window_count_list.append(window_count) # ..append window count
            window_count = 0 # ..reset window count

        if val < threshold: # If the conservation at i is less than threshold, increment count for the window
            window_count += 1

    primer_flag = False # Initialise primer flag as false

    if 1 in window_count_list[:2] or 0 in window_count_list[:2]: # If window count is at most 1 of first two windows...
        primer_flag = True # ..primer flag is true
    elif sum(window_count_list[:2]) <= 4: # If sum of window counts of the first two windows is at most 4..
        primer_flag = True # ..primer flag is true

    return primer_flag


def save_out(results, sample_id, output):
    """
    Save primer presence flags into output .txt file.

    1: primer exists
    0: primer doesn't exist
    
    First line will be the forward strand
    Second line will be the reverse strand
    """

    with open(f'{output}/{sample_id}_general_primer_out.txt', 'w') as fw:
        fw.write(f'{results[0]}\n')
        fw.write(f'{results[1]}\n')


def main():

    _PATH, _SAMPLE, _OUTPUT = parse_args()

    fwd_primer_flag = are_there_primers_in_this_sample(_PATH) # Check for general primers in fwd
    rev_primer_flag = are_there_primers_in_this_sample(_PATH, rev=True) # Check for general primers in rev

    fwd_status = '0'
    rev_status = '0'
    # Flag for primer presence: 1 for yes 0 for no
    if fwd_primer_flag:
        print('Forward primer detected!')
        fwd_status = 1
    else:
        print('No forward primer detected')
    if rev_primer_flag:
        print('Reverse primer detected!')
        rev_status = 1
    else:
        print('No reverse primer detected')

    save_out((fwd_status, rev_status), _SAMPLE, _OUTPUT) # Save primer flags to .txt file
    

if __name__ == "__main__":
    main()