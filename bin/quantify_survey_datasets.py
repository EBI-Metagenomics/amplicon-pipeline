
import os
from multiprocessing import Pool

from tqdm import tqdm

from are_there_primers import are_there_primers_in_this_sample

_PATH = '/hps/nobackup/rdf/metagenomics/service-team/users/chrisata/asv_datasets/survey_datasets/data'
_FINAL_PATH = '/hps/nobackup/rdf/metagenomics/service-team/users/chrisata/asv_datasets/survey_datasets/final_survey_dataset_data'
_CPU = 1

def assess_primer_status(samples):
    

    run_primer_count = 0

    for sample in samples:
        primer_flag = are_there_primers_in_this_sample(sample)
        if primer_flag:
            run_primer_count += 1

    if run_primer_count >= 1:
        print('Primer found!')
        return True, run_primer_count
    else:
        print('No primer found!')
        return False, run_primer_count


def main():


    total_count = 0
    primer_count = 0

    dir_list = os.listdir(_PATH)
    
    files_to_assess = []

    for dir in dir_list:
        
        if os.path.isdir(f'{_PATH}/{dir}'):
            reads_dir = f'{_PATH}/{dir}/raw'
            
            if os.path.exists(reads_dir):
                reads_list = os.listdir(reads_dir)
                
                if len(reads_list) >= 2:
                    total_count += 1
                    
                    add_list = []
                    for file in reads_list:
                        if '_1' in file or '_2' in file:
                            file_path = f'{reads_dir}/{file}'
                            add_list.append(file_path)
                    if len(add_list) == 2:
                        files_to_assess.append(add_list)
                        os.makedirs(f'{_FINAL_PATH}/{dir}/raw', exist_ok=True)
                        for file in add_list:
                            file_name = file.split('/')[-1]
                            # print(f'{_FINAL_PATH}/{dir}/raw/{file_name}')
                            os.rename(file, f'{_FINAL_PATH}/{dir}/raw/{file_name}')

                
    
    # p = Pool(_CPU) # Prepare pool for parallelisation on per sample basis    
    # out = p.map(assess_primer_status, files_to_assess)

    # primer_count = len([ res[0] for res in out if res[0] ])
    # run_primer_count = sum([ res[1] for res in out ])

    # print(files_to_assess)

    print(f'Total run count: {total_count}')
    # print(f'Total runs with primers: {primer_count}')
    # print(f'Total run files with primers: {run_primer_count}')
    # print(f'Total runs without primers: {total_count - primer_count}')



if __name__ == "__main__":
    main()