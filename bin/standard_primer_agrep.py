
import argparse
from collections import defaultdict
import os
import subprocess

from Bio.Seq import Seq
import numpy as np

from utils import primer_regex_query_builder, get_read_count

_STD_PRIMERS = "/hps/software/users/rdf/metagenomics/service-team/users/chrisata/asv_gen/data/standard_primers"

def parse_args():

    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--input", required=True, type=str, help="Path to merged FASTA to look for primers")
    parser.add_argument("-s", "--sample", required=True, type=str, help="Sample ID")
    parser.add_argument("-o", "--output", required=True, type=str, help="Output path")
    args = parser.parse_args()
  
    _INPUT = args.input
    _SAMPLE = args.sample
    _OUTPUT = args.output

    return _INPUT, _SAMPLE, _OUTPUT

def parse_std_primers():

    std_primer_dict = defaultdict(defaultdict)

    dir = os.listdir(_STD_PRIMERS)
    dir = [ f'{_STD_PRIMERS}/{path}' for path in dir ]
    
    rev_flag = False

    for path in dir:
        region = path.split('/')[-1].split('.')[0]
        with open(path, 'r') as fr:
            key = ''
            for line in fr:
                line = line.strip()
                print(region)
                if line[0] == '>':
                    if 'R' in line:
                        rev_flag = True
                    key = line[1:]
                else:
                    if rev_flag:
                        rev_conv = str(Seq(line).reverse_complement())
                        line = rev_conv
                        rev_flag = False

                    primer = primer_regex_query_builder(line)
                    std_primer_dict[region][key] = primer

    return std_primer_dict

def run_primer_agrep_once(input_path, input_primer, strand, mismatches=1):

    cmd = [
        'bash',
        '/hps/software/users/rdf/metagenomics/service-team/users/chrisata/asv_gen/bin/primer_agrep.sh',
        '-i',
        input_path,
        '-p',
        input_primer,
        '-n',
        str(mismatches),
        '-s',
        strand
    ]

    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = proc.communicate()
    output = stdout.decode('ascii')
    if output == '':
        output = 0
    else:
        output = float(output.strip())

    return output

def get_primer_props(std_primer_dict, input_path):

    threshold = 0.60
    final_primers = []

    already_searched = []
    primer_counts = []

    read_count = get_read_count(input_path, 'fastq')

    res_dict = defaultdict(defaultdict)

    for region, primer in std_primer_dict.items():
        res_dict[region]['F'] = {}
        res_dict[region]['R'] = {}

        for primer_name, primer_seq in primer.items():
            
            region_name_str = f'{region};{primer_name}'
            primer_count = 0.0

            if 'F' in primer_name:
                primer_count = run_primer_agrep_once(input_path, primer_seq, 'F', 1)
            elif 'R' in primer_name:
                primer_count = run_primer_agrep_once(input_path, primer_seq, 'R', 1)

            try:
                primer_prop = primer_count / read_count
            except ZeroDivisionError:
                primer_prop = 0

            if 'F' in primer_name:
                if primer_prop > threshold:
                    res_dict[region]['F'][primer_name] = primer_prop
            elif 'R' in primer_name:
                if primer_prop > threshold:
                    res_dict[region]['R'][primer_name] = primer_prop

            already_searched.append(region_name_str)
            print(f'{region_name_str}: {primer_prop}')
        
        if res_dict[region]['F'] == {}:
            res_dict[region].pop('F')
        if res_dict[region]['R'] == {}:
            res_dict[region].pop('R')
        

    singles = defaultdict(str)
    doubles = defaultdict(list)

    double_status = False

    for region in res_dict.keys():
        strands = res_dict[region]
                
        for strand in strands.keys():
            primers = strands[strand]
            max_prop = 0
            max_name = ''
            for primer_name, prop in primers.items():
                if prop > max_prop:
                    max_prop = prop
                    max_name = primer_name
                
            if len(strands.keys()) == 2:
                double_status = True
                doubles[region].append({max_name: max_prop})
            elif len(strands.keys()) == 1:
                singles[region] = {max_name: max_prop}

    max_region = ''
    max_primers = {}
    max_mean_prop = 0
    
    if double_status:
        for region in doubles:
            primers = doubles[region]

            f_primer_name = list(primers[0].keys())[0]
            r_primer_name = list(primers[1].keys())[0]
            f_primer_prop = primers[0][f_primer_name]
            r_primer_prop = primers[1][r_primer_name]
            
            mean_prop = (f_primer_prop + r_primer_prop) / 2.0
            if mean_prop > max_mean_prop:
                max_mean_prop = mean_prop
                max_region = region
                max_primers = [{f_primer_name: f_primer_prop}, {r_primer_name: r_primer_prop}]

    else:
        for region in singles:
            primer = singles[region]
            primer_name = list(primer.keys())[0]
            prop = primer[primer_name]
            if prop > max_mean_prop:
                max_mean_prop = prop
                max_region = region
                max_primers = {primer_name: prop}

    if max_region == '':
        print('No standard library primers!')
        return([])
    elif double_status:
        print('Standard library primers found!')
        print(f'Region: {max_region}')
        print(f'Forward Primer: {max_primers[0]}')
        print(f'Reverse Primer: {max_primers[1]}')

        return([max_region, max_primers[0], max_primers[1]])
    else:
        print('Standard library primer found on one strand!')
        print(f'Region: {max_region}')
        print(f'Primer: {max_primers}')
        
        return([max_region, max_primers])



def save_out(results, sample_id, output):

    with open(f'{output}/{sample_id}_std_primer_out.txt', 'w') as fw:
        if results == []:
            fw.write(f'')
        
        elif len(results) == 2:
            region = results[0]
            primer_name = list(results[1].keys())[0]
            primer_prop = results[1][list(results[1].keys())[0]]
            fw.write(f'{region}\n')
            fw.write(f'{primer_name}: {primer_prop}')
        
        elif len(results) == 3:
            region = results[0]
            f_primer_name = list(results[1].keys())[0]
            f_primer_prop = results[1][list(results[1].keys())[0]]
            r_primer_name = list(results[2].keys())[0]
            r_primer_prop = results[2][list(results[2].keys())[0]]
            fw.write(f'{region}\n')
            fw.write(f'{f_primer_name}: {f_primer_prop}\n')
            fw.write(f'{r_primer_name}: {r_primer_prop}')
    
def main():
    
    _INPUT, _SAMPLE, _OUTPUT = parse_args()
    std_primer_dict = parse_std_primers()
    results = get_primer_props(std_primer_dict, _INPUT)
    save_out(results, _SAMPLE, _OUTPUT)
    

if __name__ == "__main__":
    main()