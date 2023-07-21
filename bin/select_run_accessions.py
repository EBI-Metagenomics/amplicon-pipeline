
from collections import defaultdict

import pandas as pd
import numpy as np

_PATH = "/hps/nobackup/rdf/metagenomics/service-team/users/chrisata/asv_datasets/survey_datasets"

np.random.seed(12)

def main():
    
    df = pd.read_csv(f'{_PATH}/datasets.tsv', sep='\t')
    reformatted_dict = defaultdict(list)
    filtered_dict = defaultdict(list)

    for i in range(len(df)):
        run_acc = df.iloc[i, 0]
        study_acc = df.iloc[i, 1]
        
        reformatted_dict[study_acc].append(run_acc)

    for study_acc, run_list in reformatted_dict.items():

        chosen_run = np.random.choice(run_list, 1)[0]
        filtered_dict['secondary_study_accession'].append(study_acc)
        filtered_dict['run_accession'].append(chosen_run)

    filtered_df = pd.DataFrame.from_dict(filtered_dict)
    filtered_df.to_csv(f'{_PATH}/filtered_datasets.tsv', sep='\t', index=False)


if __name__ == "__main__":
    main()