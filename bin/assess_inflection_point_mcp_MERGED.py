
import argparse
from collections import defaultdict
from multiprocessing import Pool


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import tqdm

from assess_mcp_proportions_MERGED import fetch_mcp
from utils import split_dir_into_sample_paths, get_read_count, build_cons_seq

def parse_args():

    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--input", required=True, type=str, help="Path to fastq file to choose inflection point")
    parser.add_argument("-p", "--points", required=True, type=str, help="Path to inflection points file")
    parser.add_argument("-s", "--sample", required=True, type=str, help="Sample ID")
    parser.add_argument("-o", "--output", required=True, type=str, help="Output path")

    args = parser.parse_args()
    
    _PATH = args.input
    _POINTS = args.points
    _SAMPLE = args.sample
    _OUTPUT = args.output

    return _PATH, _POINTS, _SAMPLE, _OUTPUT


def assess_inflection_point_mcp_for_sample(_PATH, inf_point_list, rev=False):
    
    # TODO error handle for empty inflection point list

    beg_confs = []
    end_confs = []
    beg_cons_lens = []

    do_not_include_list = [ i + 5 for i in inf_point_list ]
    # subs_len = 15

    read_count = get_read_count(_PATH)

    n_prop = 0.95

    for beg in inf_point_list:
        beg += 4
        mcp_count_dict = fetch_mcp(_PATH, beg, rev=rev)
        mcp_cons_list = []
        mcp_len = len(list(mcp_count_dict.keys())[0])

        for i in range(mcp_len):
            index_base_dict = defaultdict(int)
            for mcp in mcp_count_dict.keys():
                if len(mcp) < mcp_len:
                    continue
                base = mcp[i]
                index_base_dict[base] += mcp_count_dict[mcp]
            mcp_cons_list.append(index_base_dict)

        cons_seq, cons_confs = build_cons_seq(mcp_cons_list, read_count, n_prop, do_not_include_list)
        mcp_sum = sum(mcp_count_dict.values())
        prop = mcp_sum/read_count

        n_count = float(cons_seq.count('N'))/len(cons_seq)

        beg_confs.append(np.mean(cons_confs))
        beg_cons_lens.append(len(cons_seq))
        # fwd_beg_confs.append(n_count)
        # fwd_beg_confs.append(fwd_prop)


        # fwd_confs.append(fwd_prop)
        print(1, beg, cons_seq)

    for i, end in enumerate(inf_point_list):
        end += 5
        subs_len = beg_cons_lens[i]
        l = end + subs_len - 1
        mcp_count_dict = fetch_mcp(_PATH, l, end)
        mcp_cons_list = []
        mcp_len = len(list(mcp_count_dict.keys())[0])

        for i in range(mcp_len):
            index_base_dict = defaultdict(int)
            for mcp in mcp_count_dict.keys():
                if len(mcp) < subs_len or len(mcp) < mcp_len:
                    continue
                base = mcp[i]
                index_base_dict[base] += mcp_count_dict[mcp]
            mcp_cons_list.append(index_base_dict)

        cons_seq, cons_confs = build_cons_seq(mcp_cons_list, read_count, n_prop, do_not_include_list, end)
        mcp_sum = sum(mcp_count_dict.values())
        prop = mcp_sum/read_count

        n_count = float(cons_seq.count('N'))/len(cons_seq)

        end_confs.append(np.mean(cons_confs))
        # fwd_end_confs.append(n_count)
        # fwd_end_confs.append(fwd_prop)

        # fwd_confs.append(fwd_prop)
        print(end, l, cons_seq)

    diff_res = [ beg_confs[i] - end_confs[i] for i in range(len(beg_confs))]
    diff_res_sorted = sorted(diff_res, reverse=True)

    print(beg_confs)
    print(end_confs)
    print(diff_res)

    ini_max_res = diff_res_sorted[0]
    curr_max_index = diff_res.index(ini_max_res)

    for res in diff_res_sorted[1:]:
        curr_res_index = np.where(diff_res == res)[0][0]

        index_diff = inf_point_list[curr_max_index] - inf_point_list[curr_res_index]

        if ini_max_res - res < 0.05 and ( index_diff <= 3 and index_diff > 0 ):
            curr_max_index = curr_res_index

    cutoff = inf_point_list[curr_max_index] + 5

    print(cutoff)

    return cutoff

def main():

    _PATH, _POINTS, _SAMPLE, _OUTPUT = parse_args()
    inf_df = pd.read_csv(_POINTS, sep='\t')

    f_slice = inf_df[inf_df.strand == 'F']
    r_slice = inf_df[inf_df.strand == 'R']
    r_slice = r_slice.reset_index(drop=True)

    if not f_slice.empty:
        inf_list = f_slice.inf_point.tolist()
        print(inf_list)
        assess_inflection_point_mcp_for_sample(_PATH, inf_list)

    if not r_slice.empty:
        inf_list = r_slice.inf_point.tolist()
        print(inf_list)
        assess_inflection_point_mcp_for_sample(_PATH, inf_list, rev=True)


if __name__ == "__main__":
    main()