
import argparse
from collections import defaultdict

import numpy as np

from assess_mcp_proportions_MERGED import fetch_mcp
from utils import get_read_count, build_cons_seq

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

    read_count = get_read_count(_PATH, 'fastq')
    subs_len = 100

    mcp_count_dict = fetch_mcp(_PATH, subs_len, rev=rev)
    mcp_cons_list = []

    for i in range(subs_len):
        index_base_dict = defaultdict(int)
        for mcp in mcp_count_dict.keys():
            if len(mcp) < subs_len:
                continue
            base = mcp[i]
            index_base_dict[base] += mcp_count_dict[mcp]
        mcp_cons_list.append(index_base_dict)

    cons_seq, cons_confs = build_cons_seq(mcp_cons_list, read_count, cons_threshold=0.9)

    drops_list = []
    drops_list_post_cut = []
    final_window = []
    window_size = 10
    window_count = 0
    window_count_post_cut = 0
    test_window = []
    test_window_post_cut = []

    max_cons = np.quantile(cons_confs, 0.75)
    threshold = max_cons - 0.15
    if max_cons < 0.75:
        threshold = 0.75
    max_cons_post_cut = np.quantile(cons_confs[20:], 0.75)
    threshold_post_cut = max_cons_post_cut - 0.15
    if max_cons_post_cut < 0.75:
        threshold_post_cut = 0.75

    if max_cons < 0.6:
        return False

    for i, val in enumerate(cons_confs):
        if i%window_size == 0 and i !=0:
            subwindow = [window_count] * window_size
            drops_list.extend(subwindow)
            final_window.append(window_count)
            window_count = 0
            if i > 20:
                subwindow = [window_count_post_cut] * window_size
                drops_list_post_cut.extend(subwindow)
                window_count_post_cut = 0 
        if val < threshold:
            window_count += 1
            test_window.append(1)
            if i > 20:
                test_window_post_cut.append(1)
        else:
            test_window.append(0)
            if i > 20:
                test_window_post_cut.append(0)
        if i > 20 and val < threshold_post_cut:
            window_count_post_cut += 1
        

    subwindow = [window_count] * window_size
    drops_list.extend(subwindow)
    subwindow = [window_count_post_cut] * window_size
    drops_list_post_cut.extend(subwindow)

    primer_flag = False

    if 1 in final_window[:2] or 0 in final_window[:2]:
        primer_flag = True
    elif sum(final_window[:2]) <= 4:
        primer_flag = True

    return primer_flag


def save_out(results, sample_id, output):

    with open(f'{output}/{sample_id}_general_primer_out.txt', 'w') as fw:
        fw.write(f'{results[0]}\n')
        fw.write(f'{results[1]}\n')


def main():

    _PATH, _SAMPLE, _OUTPUT = parse_args()

    fwd_primer_flag = are_there_primers_in_this_sample(_PATH)
    rev_primer_flag = are_there_primers_in_this_sample(_PATH, True)

    fwd_status = '0'
    rev_status = '0'

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

    save_out((fwd_status, rev_status), _SAMPLE, _OUTPUT)
    

if __name__ == "__main__":
    main()