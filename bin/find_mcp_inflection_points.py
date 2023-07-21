
from collections import defaultdict
import pickle
import argparse


import pandas as pd
import numpy as np
from scipy import signal
from scipy.ndimage import gaussian_filter1d


import matplotlib.pyplot as plt
import seaborn as sns


# fwd_path = '/hps/nobackup/rdf/metagenomics/service-team/users/chrisata/asv_datasets/PRJNA714811/raw/mcp_out/fwd_mcp_props_w4.pkl'
# rev_path = '/hps/nobackup/rdf/metagenomics/service-team/users/chrisata/asv_datasets/PRJNA714811/raw/mcp_out/rev_mcp_props_w4.pkl'
# fwd_path = '/hps/nobackup/rdf/metagenomics/service-team/users/chrisata/asv_datasets/PRJNA714811/raw/mcp_out/fwd_mcp_cons_w4.pkl'
# rev_path = '/hps/nobackup/rdf/metagenomics/service-team/users/chrisata/asv_datasets/PRJNA714811/raw/mcp_out/rev_mcp_cons_w4.pkl'

# fwd_path = '/hps/nobackup/rdf/metagenomics/service-team/users/chrisata/asv_datasets/SRP386045/raw/mcp_out/fwd_mcp_props_w4.pkl'
# rev_path = '/hps/nobackup/rdf/metagenomics/service-team/users/chrisata/asv_datasets/SRP386045/raw/mcp_out/rev_mcp_props_w4.pkl'
# fwd_path = '/hps/nobackup/rdf/metagenomics/service-team/users/chrisata/asv_datasets/SRP386045/raw/mcp_out/fwd_mcp_cons_w4.pkl'
# rev_path = '/hps/nobackup/rdf/metagenomics/service-team/users/chrisata/asv_datasets/SRP386045/raw/mcp_out/rev_mcp_cons_w4.pkl'


# fwd_path = '/hps/nobackup/rdf/metagenomics/service-team/users/chrisata/asv_datasets/ERP123542/raw/mcp_out/fwd_mcp_props_w4.pkl'
# rev_path = '/hps/nobackup/rdf/metagenomics/service-team/users/chrisata/asv_datasets/ERP123542/raw/mcp_out/rev_mcp_props_w4.pkl'
# fwd_path = '/hps/nobackup/rdf/metagenomics/service-team/users/chrisata/asv_datasets/ERP123542/raw/mcp_out/fwd_mcp_cons_w4.pkl'
# rev_path = '/hps/nobackup/rdf/metagenomics/service-team/users/chrisata/asv_datasets/ERP123542/raw/mcp_out/rev_mcp_cons_w4.pkl'

# fwd_path = '/hps/nobackup/rdf/metagenomics/service-team/users/chrisata/asv_datasets/SRP293247/raw/mcp_out/fwd_mcp_props_w4.pkl'
# rev_path = '/hps/nobackup/rdf/metagenomics/service-team/users/chrisata/asv_datasets/SRP293247/raw/mcp_out/rev_mcp_props_w4.pkl'
# fwd_path = '/hps/nobackup/rdf/metagenomics/service-team/users/chrisata/asv_datasets/SRP293247/raw/mcp_out/fwd_mcp_cons_w4.pkl'
# rev_path = '/hps/nobackup/rdf/metagenomics/service-team/users/chrisata/asv_datasets/SRP293247/raw/mcp_out/rev_mcp_cons_w4.pkl'

# fwd_path = '/hps/nobackup/rdf/metagenomics/service-team/users/chrisata/asv_datasets/ERP124393/raw/mcp_out/fwd_mcp_props_w4.pkl'
# rev_path = '/hps/nobackup/rdf/metagenomics/service-team/users/chrisata/asv_datasets/ERP124393/raw/mcp_out/rev_mcp_props_w4.pkl'
# fwd_path = '/hps/nobackup/rdf/metagenomics/service-team/users/chrisata/asv_datasets/ERP124393/raw/mcp_out/fwd_mcp_cons_w4.pkl'
# rev_path = '/hps/nobackup/rdf/metagenomics/service-team/users/chrisata/asv_datasets/ERP124393/raw/mcp_out/rev_mcp_cons_w4.pkl'


parser = argparse.ArgumentParser()

parser.add_argument("-d", "--dir", required=True, type=str, help="Input directory containing mcp pickle files")
parser.add_argument("-cpu", "--cpu", type=int, default=1, help="Number of CPUs")

args = parser.parse_args()
_DIR = args.dir
_CPU = args.cpu


def main():

    # fwd_path = f'{_DIR}/fwd_mcp_props_w4.pkl'
    # rev_path = f'{_DIR}/rev_mcp_props_w4.pkl'
    fwd_path = f'{_DIR}/fwd_mcp_cons_w4.pkl'
    rev_path = f'{_DIR}/rev_mcp_cons_w4.pkl'
    out_path = f'{_DIR}/inf_points'

    with open(fwd_path, 'rb') as f:
        fwd_df = pickle.load(f)

    with open(rev_path, 'rb') as f:
        rev_df = pickle.load(f)

    x = [ int(i) for i in fwd_df.columns.tolist() ]

    fwd_len_dict = defaultdict(int)
    rev_len_dict = defaultdict(int)
    subs_len = 4


    for i in range(len(fwd_df)):
        arr_fwd = fwd_df.iloc[i].tolist()
        arr_fwd = [ -val for val in arr_fwd ]
        sample = fwd_df.index.tolist()[i]


        m_fwd = np.diff(arr_fwd)/np.diff(x)
        n_fwd = np.diff(m_fwd)

        infl_points_fwd = np.where(m_fwd > np.percentile(m_fwd, 80))[0]
        # infl_points_fwd = np.argsort(m_fwd)[-5:]


        arr_rev = rev_df.iloc[i].tolist()
        arr_rev = [ -val for val in arr_rev ]

        m_rev = np.diff(arr_rev)/np.diff(x)
        n_rev = np.diff(m_rev)
        infl_points_rev = np.where(m_rev > np.percentile(m_rev, 80))[0]
        # infl_points_rev = np.argsort(m_rev)[-5:]

        with open(f'{out_path}/{sample}_fwd_inf_points.txt', 'w') as f:
            for ind in infl_points_fwd:
                inf_point = x[ind]
                # if inf_point < 10:
                #     continue
                fwd_len_dict[inf_point] += 1
                f.write(f'{inf_point}\n')

                
        with open(f'{out_path}/{sample}_rev_inf_points.txt', 'w') as f:
            for ind in infl_points_rev:
                inf_point = x[ind]
                # if inf_point < 10:
                #     continue
                rev_len_dict[inf_point] += 1
                f.write(f'{inf_point}\n')

    print(x)
    print(fwd_len_dict)
    print(rev_len_dict)


if __name__ == "__main__":
    main()