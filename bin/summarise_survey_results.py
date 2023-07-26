
from collections import defaultdict
import os

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

_RES_PATH = "/hps/nobackup/rdf/metagenomics/service-team/users/chrisata/asv_datasets/survey_datasets/final_survey_dataset_data/merged/merged"
_OUT_PATH = "/hps/nobackup/rdf/metagenomics/service-team/users/chrisata/asv_datasets/survey_datasets/final_survey_dataset_data/summary"

def parse_general_out(path):

    try:
        with open(path, 'r') as fr:
            lines = list(fr)
            res = [ line.strip() for line in lines ]
    except FileNotFoundError:
        return 0
    
    return res

def parse_std_out(path):

    try:
        with open(path, 'r') as fr:
            lines = list(fr)

            if lines == []:
                return '', '', ''
            
            res = [ line.strip() for line in lines ]
            primer_region = res[0]
            
            if len(res) == 2:
                
                primer_name = res[1].split(':')[0]

                if 'F' in res[1]:
                    return primer_region, primer_name, ''
                elif 'R' in res[1]:
                    return primer_region, '', primer_name
            
            elif len(res) == 3:
                
                fwd_primer_name = res[1].split(':')[0]
                rev_primer_name = res[2].split(':')[0]

                return primer_region, fwd_primer_name, rev_primer_name        

    except FileNotFoundError:
        return 0, 0, 0
    
    return res

def main():
    
    projs = os.listdir(_RES_PATH)

    res_dict = defaultdict(list)

    for proj in projs:
        
        proj_path = f'{_RES_PATH}/{proj}'
        proj_files = os.listdir(proj_path)
        run_id = proj_files[0].split('_')[0]

        general_out_path = f'{proj_path}/{run_id}_MERGED_general_primer_out.txt'
        std_out_path = f'{proj_path}/{run_id}_MERGED_std_primer_out.txt'
        general_res = parse_general_out(general_out_path)
        std_region, std_fwd, std_rev = parse_std_out(std_out_path)

        if general_res == 0 or std_region == 0:
            continue
        
        std_res = []
        if std_fwd == '':
            std_res.append('0')
        else:
            std_res.append('1')
        if std_rev == '':
            std_res.append('0')
        else:
            std_res.append('1')
        
        res_dict['proj'].append(proj)
        res_dict['general_fwd'].append(general_res[0])
        res_dict['general_rev'].append(general_res[1])
        res_dict['general_comb'].append('-'.join(general_res))
        res_dict['std_fwd'].append(std_res[0])
        res_dict['std_rev'].append(std_res[1])
        res_dict['std_comb'].append('-'.join(std_res))
        if std_fwd != '':
            res_dict['std_fwd_primer_region'].append(std_region)
            res_dict['std_fwd_primer_name'].append(std_fwd)
        else:
            res_dict['std_fwd_primer_region'].append('none')
            res_dict['std_fwd_primer_name'].append('none')
        if std_rev != '':
            res_dict['std_rev_primer_region'].append(std_region)
            res_dict['std_rev_primer_name'].append(std_rev)
        else:
            res_dict['std_rev_primer_region'].append('none')
            res_dict['std_rev_primer_name'].append('none')
        
    
    res_df = pd.DataFrame.from_dict(res_dict)

    print(res_df)

    res_counts = defaultdict(int)

    for i in range(len(res_df)):

        general_res = res_df.iloc[i, 3]
        std_res = res_df.iloc[i, 6]

        if general_res == '1-1' and std_res == '1-1':
            res_counts['GGSS'] += 1
        elif general_res == '0-0' and std_res == '0-0':
            res_counts['ggss'] += 1
        elif general_res == '1-0' and std_res == '1-0':
            res_counts['GgSs'] += 1
        elif general_res == '0-1' and std_res == '0-1':
            res_counts['gGsS'] += 1
        elif general_res == '1-1' and (std_res == '0-1' or std_res == '1-0' or std_res == '0-0'):
            # print(res_df.iloc[i, 0], general_res, std_res)
            res_counts['GG__'] += 1
        elif general_res == '1-0' and std_res == '0-0':
            # print(res_df.iloc[i, 0], general_res, std_res)
            res_counts['GG__'] += 1
        elif general_res == '0-1' and std_res == '0-0':
            # print(res_df.iloc[i, 0], general_res, std_res)
            res_counts['GG__'] += 1
        elif general_res == '0-0' and (std_res == '0-1' or std_res == '1-0' or std_res == '1-1'):
            # print(res_df.iloc[i, 0], general_res, std_res)
            res_counts['gg__'] += 1
        elif general_res == '1-0' and (std_res == '0-1' or std_res == '1-1'):
            # print(res_df.iloc[i, 0], general_res, std_res)
            res_counts['gg__'] += 1
        elif general_res == '0-1' and (std_res == '1-0' or std_res == '1-1'):
            # print(res_df.iloc[i, 0], general_res, std_res)
            res_counts['gg__'] += 1


    res_counts_df = pd.DataFrame(res_counts, index=[0])
    print(res_counts_df)

    filtered_res_df = res_df[(res_df.general_comb == '1-1') & (res_df.std_comb == '1-1')]
    long_res_counts_df = pd.melt(res_counts_df)

    sns.set(font_scale=0.25)
    
    sns.barplot(long_res_counts_df, x='variable', y='value')
    plt.savefig(f'{_OUT_PATH}/flag_counts_barplot.png', dpi=200)

    plt.figure() # reset plot

    sns.countplot(data=filtered_res_df, x='std_fwd_primer_region')
    plt.savefig(f'{_OUT_PATH}/std_fwd_primer_region_countplot.png', dpi=200)

    plt.figure() # reset plot

    sns.countplot(data=filtered_res_df, x='std_fwd_primer_name')
    plt.savefig(f'{_OUT_PATH}/std_fwd_primer_name_countplot.png', dpi=200)

    plt.figure() # reset plot

    sns.countplot(data=filtered_res_df, x='std_rev_primer_name')
    plt.savefig(f'{_OUT_PATH}/std_rev_primer_name_countplot.png', dpi=200)


if __name__ == "__main__":
    main()