
from collections import defaultdict
import argparse

import pandas as pd
import numpy as np

def parse_args():

    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--input", required=True, type=str, help="Path to mcp tsv file to find inflection points")
    parser.add_argument("-s", "--sample", required=True, type=str, help="Sample ID")
    parser.add_argument("-o", "--output", required=True, type=str, help="Output path")

    args = parser.parse_args()
    
    _PATH = args.input
    _SAMPLE = args.sample
    _OUTPUT = args.output

    return _PATH, _SAMPLE, _OUTPUT

def main():

    _PATH, _SAMPLE, _OUTPUT = parse_args()

    mcp_df = pd.read_csv(_PATH, sep='\t', index_col=0)
    inf_point_dict = defaultdict(list)

    x = [ int(i) for i in mcp_df.columns.tolist() ]

    for i in range(len(mcp_df)):
        strand = mcp_df.index[i]
        props = mcp_df.iloc[i].tolist()
        props = [ -val for val in props ]

        prop_diff = np.diff(props)/np.diff(x)

        infl_points = np.where(prop_diff > np.percentile(prop_diff, 80))[0]
        for ind in infl_points:
            inf_point = x[ind]

            if inf_point < 15 or inf_point > 20:
                continue

            inf_point_dict['strand'].append(strand)
            inf_point_dict['inf_point'].append(inf_point)

    if len(inf_point_dict) > 0:
        inf_point_df = pd.DataFrame.from_dict(inf_point_dict)
        inf_point_df.to_csv(f'{_OUTPUT}/{_SAMPLE}_inf_points.tsv', sep='\t', index=False)


if __name__ == "__main__":
    main()