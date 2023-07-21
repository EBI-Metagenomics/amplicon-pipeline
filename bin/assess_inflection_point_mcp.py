
import argparse
from collections import defaultdict
from multiprocessing import Pool


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import tqdm

from assess_mcp_proportions import fetch_mcp
from utils import split_dir_into_sample_paths, get_read_count, build_cons_seq


parser = argparse.ArgumentParser()

parser.add_argument("-d", "--dir", required=True, type=str, help="Input directory with both strands for every sample")
parser.add_argument("-cpu", "--cpu", type=int, default=1, help="Number of CPUs")

args = parser.parse_args()
_DIR = args.dir
_CPU = args.cpu

def extract_inf_points(fileroot):
    
    inf_dir_path = f'{_DIR}/mcp_out/inf_points'
    fwd_inf_point_list = []
    rev_inf_point_list = []


    with open(f'{inf_dir_path}/{fileroot}_fwd_inf_points.txt', 'r') as f:
        fwd_inf_point_list = f.readlines()
        fwd_inf_point_list = sorted([ int(line.strip()) for line in fwd_inf_point_list ])

    with open(f'{inf_dir_path}/{fileroot}_rev_inf_points.txt', 'r') as f:
        rev_inf_point_list = f.readlines()
        rev_inf_point_list = sorted([ int(line.strip()) for line in rev_inf_point_list ])

    return fwd_inf_point_list, rev_inf_point_list


def assess_inflection_point_mcp_for_sample(sample):
    
    fileroot = sample.split('/')[-1]
    fwd_inf_point_list, rev_inf_point_list = extract_inf_points(fileroot)

    if len(fwd_inf_point_list) == 0 or len(rev_inf_point_list) == 0:
        return 0, 0


    # fwd_inf_point_list = [15, 16, 17, 18, 19, 20, 23]
    # rev_inf_point_list = [18, 19, 22, 23]


    # fwd_inf_point_list = list(range(2, 25, 1))
    # rev_inf_point_list = list(range(2, 25, 1))
    
    print(f'Processing {fileroot}')
    # TODO: need to automate fileroots 
    fwd = f"{sample}_1.fastq.gz"
    rev = f"{sample}_2.fastq.gz"

    fwd_beg_confs = []
    rev_beg_confs = []
    fwd_end_confs = []
    rev_end_confs = []
    fwd_beg_cons_lens = []
    rev_beg_cons_lens = []

    fwd_do_not_include = [ i + 5 for i in fwd_inf_point_list ]
    rev_do_not_include = [ i + 5 for i in rev_inf_point_list ]
    # subs_len = 15

    # fwd_inf_point_list = [ i for i in fwd_inf_point_list if i >=8]
    # rev_inf_point_list = [ i for i in rev_inf_point_list if i >=8]

    fwd_read_count = get_read_count(fwd)
    rev_read_count = get_read_count(rev)


    n_prop = 0.95

    # print('Forward strand')
    for end in fwd_inf_point_list:
        end += 4
        fwd_mcp_count_dict = fetch_mcp(fwd, end)
        fwd_mcp_cons_list = []
        mcp_len = len(list(fwd_mcp_count_dict.keys())[0])

        for i in range(mcp_len):
            index_base_dict = defaultdict(int)
            for mcp in fwd_mcp_count_dict.keys():
                if len(mcp) < mcp_len:
                    continue
                base = mcp[i]
                index_base_dict[base] += fwd_mcp_count_dict[mcp]
            fwd_mcp_cons_list.append(index_base_dict)

        cons_seq, cons_confs = build_cons_seq(fwd_mcp_cons_list, fwd_read_count, n_prop, fwd_do_not_include)
        fwd_mcp_sum = sum(fwd_mcp_count_dict.values())
        fwd_prop = fwd_mcp_sum/fwd_read_count

        n_count = float(cons_seq.count('N'))/len(cons_seq)

        fwd_beg_confs.append(cons_confs)
        fwd_beg_cons_lens.append(len(cons_seq))
        # fwd_beg_confs.append(n_count)
        # fwd_beg_confs.append(fwd_prop)


        # fwd_confs.append(fwd_prop)
        print(1, end, cons_seq)

    # print('Reverse strand')
    for end in rev_inf_point_list:
        end += 4
        rev_mcp_count_dict = fetch_mcp(rev, end)

        rev_mcp_cons_list = []
        mcp_len = len(list(rev_mcp_count_dict.keys())[0])

        for i in range(mcp_len):
            index_base_dict = defaultdict(int)
            for mcp in rev_mcp_count_dict.keys():
                if len(mcp) < mcp_len:
                    continue
                base = mcp[i]
                index_base_dict[base] += rev_mcp_count_dict[mcp]
            rev_mcp_cons_list.append(index_base_dict)

        cons_seq, cons_confs = build_cons_seq(rev_mcp_cons_list, rev_read_count, n_prop, rev_do_not_include)
        rev_mcp_sum = sum(rev_mcp_count_dict.values())
        rev_prop = rev_mcp_sum/rev_read_count

        n_count = float(cons_seq.count('N'))/len(cons_seq)

        rev_beg_confs.append(cons_confs)
        rev_beg_cons_lens.append(len(cons_seq))
        # rev_beg_confs.append(n_count)
        # rev_beg_confs.append(rev_prop)

        # rev_confs.append(rev_prop)
        print(1, end, cons_seq)



    # print('Forward strand')
    for i, beg in enumerate(fwd_inf_point_list):
        beg += 5
        subs_len = fwd_beg_cons_lens[i]
        l = beg + subs_len
        fwd_mcp_count_dict = fetch_mcp(fwd, l, beg)
        fwd_mcp_cons_list = []
        mcp_len = len(list(fwd_mcp_count_dict.keys())[0])

        for i in range(mcp_len):
            index_base_dict = defaultdict(int)
            for mcp in fwd_mcp_count_dict.keys():
                if len(mcp) < subs_len or len(mcp) < mcp_len:
                    continue
                base = mcp[i]
                index_base_dict[base] += fwd_mcp_count_dict[mcp]
            fwd_mcp_cons_list.append(index_base_dict)

        cons_seq, cons_confs = build_cons_seq(fwd_mcp_cons_list, fwd_read_count, n_prop, fwd_do_not_include, beg)
        fwd_mcp_sum = sum(fwd_mcp_count_dict.values())
        fwd_prop = fwd_mcp_sum/fwd_read_count

        n_count = float(cons_seq.count('N'))/len(cons_seq)

        fwd_end_confs.append(cons_confs)
        # fwd_end_confs.append(n_count)
        # fwd_end_confs.append(fwd_prop)

        # fwd_confs.append(fwd_prop)
        print(beg, l, cons_seq)

    # print('Reverse strand')
    for i, beg in enumerate(rev_inf_point_list):
        beg += 5
        subs_len = rev_beg_cons_lens[i]
        l = beg + subs_len
        rev_mcp_count_dict = fetch_mcp(rev, l, beg)

        rev_mcp_cons_list = []
        mcp_len = len(list(rev_mcp_count_dict.keys())[0])

        for i in range(mcp_len):
            index_base_dict = defaultdict(int)
            for mcp in rev_mcp_count_dict.keys():
                if len(mcp) < subs_len or len(mcp) < mcp_len:
                    continue
                base = mcp[i]
                index_base_dict[base] += rev_mcp_count_dict[mcp]
            rev_mcp_cons_list.append(index_base_dict)

        cons_seq, cons_confs = build_cons_seq(rev_mcp_cons_list, rev_read_count, n_prop, rev_do_not_include, beg)
        rev_mcp_sum = sum(rev_mcp_count_dict.values())
        rev_prop = rev_mcp_sum/rev_read_count

        n_count = float(cons_seq.count('N'))/len(cons_seq)


        rev_end_confs.append(cons_confs)
        # rev_end_confs.append(n_count)
        # rev_end_confs.append(rev_prop)
        # rev_confs.append(rev_prop)
        print(beg, l, cons_seq)





    # fwd_res = [ fwd_end_confs[i] - fwd_beg_confs[i] for i in range(len(fwd_beg_confs))]
    # fwd_res_sorted = sorted(fwd_res, reverse=True)
    # rev_res = [ rev_end_confs[i] - rev_beg_confs[i] for i in range(len(rev_beg_confs))]
    # rev_res_sorted = sorted(rev_res, reverse=True)


    fwd_res = [ fwd_beg_confs[i] - fwd_end_confs[i] for i in range(len(fwd_beg_confs))]
    fwd_res_sorted = sorted(fwd_res, reverse=True)
    rev_res = [ rev_beg_confs[i] - rev_end_confs[i] for i in range(len(rev_beg_confs))]
    rev_res_sorted = sorted(rev_res, reverse=True)


    print(fwd_beg_confs)
    print(fwd_end_confs)
    print(fwd_res)
    print(rev_beg_confs)
    print(rev_end_confs)
    print(rev_res)

    m_fwd = np.diff(fwd_res)/np.diff(list(range(len(fwd_inf_point_list))))
    m_rev = np.diff(rev_res)/np.diff(list(range(len(rev_inf_point_list))))

    # fwd_max = np.max(m_fwd)
    # rev_max = np.max(m_rev)
    # fwd_ind = np.where(m_fwd == fwd_max)[0][0]
    # rev_ind = np.where(m_rev == rev_max)[0][0]
    # fwd_curr_max_index = fwd_ind
    # rev_curr_max_index = rev_ind

    # fwd_res_sorted = sorted(m_fwd, reverse=True)
    # rev_res_sorted = sorted(m_rev, reverse=True)

    # print(m_fwd)
    # print(m_rev)
    # print(fwd_ind)
    # print(rev_ind)

    fwd_diffs = []

    for i in range(len(fwd_res) - 1):
        beg = fwd_res[i]
        end = fwd_res[i+1]

        diff = end - beg
        fwd_diffs.append(diff)


    rev_diffs = []

    for i in range(len(rev_res) - 1):
        beg = rev_res[i]
        end = rev_res[i+1]

        diff = end - beg
        rev_diffs.append(diff)


    fwd_max_val = np.max(fwd_diffs)
    rev_max_val = np.max(rev_diffs)
    fwd_max_index = np.where(fwd_diffs == np.max(fwd_diffs))[0][0]
    rev_max_index = np.where(rev_diffs == np.max(rev_diffs))[0][0]

    fwd_diffs_sorted = sorted(fwd_diffs)
    rev_diffs_sorted = sorted(rev_diffs)

    for res in fwd_diffs_sorted[1:]:
        curr_res_index = np.where(fwd_diffs == res)[0][0]
        # if res < 0:
        #     continue
        # and curr_res_index < fwd_curr_max_index:
        max_potential_cutoff = fwd_inf_point_list[fwd_max_index + 1] + 5
        curr_potential_cutoff = fwd_inf_point_list[curr_res_index + 1] + 5
        if fwd_max_val - res < 0.03 and (((max_potential_cutoff - curr_potential_cutoff) <= 3) and ((max_potential_cutoff - curr_potential_cutoff) > 0)):
            # fwd_curr_max_res = res
            fwd_max_index = curr_res_index

    for res in rev_diffs_sorted[1:]:
        curr_res_index = np.where(rev_diffs == res)[0][0]
        # if res < 0:
        #     continue
        # and curr_res_index < fwd_curr_max_index:
        max_potential_cutoff = rev_inf_point_list[rev_max_index + 1] + 5
        curr_potential_cutoff = rev_inf_point_list[curr_res_index + 1] + 5
        if rev_max_val - res < 0.03 and (((max_potential_cutoff - curr_potential_cutoff) <= 3) and ((max_potential_cutoff - curr_potential_cutoff) > 0)):
            # fwd_curr_max_res = res
            rev_max_index = curr_res_index


    print(fwd_diffs)
    print(rev_diffs)

  
    # ini_fwd_curr_max_res = fwd_res_sorted[0]
    # # # fwd_curr_max_res = fwd_res_sorted[0]
    # fwd_curr_max_index = fwd_res.index(ini_fwd_curr_max_res)

    # for res in fwd_res_sorted[1:]:
    #     curr_res_index = np.where(fwd_res == res)[0][0]
    #     # if res < 0:
    #     #     continue
    #     # and curr_res_index < fwd_curr_max_index:
    #     if ini_fwd_curr_max_res - res < 0.05 and curr_res_index < fwd_curr_max_index:
    #         # fwd_curr_max_res = res
    #         fwd_curr_max_index = curr_res_index

    # ini_rev_curr_max_res = rev_res_sorted[0]
    # # # fwd_curr_max_res = fwd_res_sorted[0]
    # rev_curr_max_index = rev_res.index(ini_rev_curr_max_res)


    # for res in rev_res_sorted[1:]:
    #     curr_res_index = np.where(rev_res == res)[0][0]
    #     # if res < 0:
    #     #     continue
    #     if ini_rev_curr_max_res - res < 0.05 and curr_res_index < rev_curr_max_index:
    #         # fwd_curr_max_res = res
    #         rev_curr_max_index = curr_res_index


    # ini_rev_curr_max_res = rev_res_sorted[0]
    # # rev_curr_max_res = rev_res_sorted[0]
    # rev_curr_max_index = rev_res.index(ini_rev_curr_max_res)

    # for res in rev_res_sorted[1:]:
    #     curr_res_index = rev_res.index(res)
    #     # if res < 0:
    #     #     continue
    #     if ini_rev_curr_max_res - res < 0.03 and curr_res_index < rev_curr_max_index:
    #         # rev_curr_max_res = res
    #         rev_curr_max_index = curr_res_index


    # print(fwd_end_confs)
    # print(fwd_beg_confs)

    # print(rev_end_confs)
    # print(rev_beg_confs)

    # print(fwd_res)
    # print(rev_res)

    # print(fwd_curr_max_index)
    # print(rev_curr_max_index)




    fwd_cutoff = fwd_inf_point_list[fwd_max_index + 1] + 5
    rev_cutoff = rev_inf_point_list[rev_max_index + 1] + 5

    print(fwd_cutoff, rev_cutoff)

    return fwd_cutoff, rev_cutoff

    # plt.plot(fwd_inf_point_list, fwd_confs)
    # plt.ylim(0, 1)
    # plt.axvline(x=5, color='red')
    # plt.axvline(x=9, color='red')
    # plt.axvline(x=15, color='red')
    # plt.axvline(x=19, color='green')
    # plt.savefig(f"/hps/nobackup/rdf/metagenomics/service-team/users/chrisata/asv_datasets/PRJNA714811/raw/mcp_out/PRJNA714811_test_{subs_len}_inflection_fwd.png", dpi=200)

    # plt.figure()

    # plt.plot(rev_inf_point_list, rev_confs)
    # plt.ylim(0, 1)
    # plt.axvline(x=3, color='red')
    # plt.axvline(x=4, color='red')
    # plt.axvline(x=17, color='red')
    # plt.axvline(x=21, color='green')
    # plt.savefig(f"/hps/nobackup/rdf/metagenomics/service-team/users/chrisata/asv_datasets/PRJNA714811/raw/mcp_out/PRJNA714811_test_{subs_len}_inflection_rev.png", dpi=200)


def main():

    sample_list = split_dir_into_sample_paths(_DIR)
    sample_names = [ sample.split('/')[-1] for sample in sample_list ]

    res_dict = defaultdict(list)


    p = Pool(_CPU) # Prepare pool for parallelisation on per sample basis

    out = p.map(assess_inflection_point_mcp_for_sample, sample_list)

    for res in out:
        fwd_cutoff = res[0]
        rev_cutoff = res[1]
        # fwd_cutoff, rev_cutoff = assess_inflection_point_mcp_for_sample(sample)

        res_dict['inf_cutoff'].append(fwd_cutoff)
        res_dict['strand'].append('fwd')
        res_dict['inf_cutoff'].append(rev_cutoff)
        res_dict['strand'].append('rev')
        

    res_df = pd.DataFrame.from_dict(res_dict)
    print(res_df)

    sns.stripplot(data=res_df, x='inf_cutoff', y='strand', hue='strand', alpha=0.75)
    plt.axvline(x=18, ymin=0.5, color='red')
    plt.axvline(x=22, ymax=0.5, color='green')
    plt.xlim(7, 30)
    plt.xticks(range(7, 30))
    
    plt.savefig(f'{_DIR}/mcp_out/inf_cutoff_selection_plot_cons_test.png', dpi=200)


if __name__ == "__main__":
    main()