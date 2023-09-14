
import argparse
from collections import defaultdict
import os
import subprocess

from Bio.Seq import Seq

from utils import primer_regex_query_builder, get_read_count

# Folder containing the library of standard primers as fasta files
_STD_PRIMERS = "/hps/software/users/rdf/metagenomics/service-team/users/chrisata/asv_gen/data/standard_primers"

def parse_args():

    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--input", required=True, type=str, help="Path to merged FASTQ to look for primers")
    parser.add_argument("-s", "--sample", required=True, type=str, help="Sample ID")
    parser.add_argument("-o", "--output", required=True, type=str, help="Output path")
    args = parser.parse_args()
  
    _INPUT = args.input
    _SAMPLE = args.sample
    _OUTPUT = args.output

    return _INPUT, _SAMPLE, _OUTPUT

def parse_std_primers():
    """
    Parse the library of standard primers.

    Reads the fasta files in the given directory "_STD_PRIMERS"
    Primer names (which are the fasta headers) are labeled with F or R for 5'-3' and 3'-5' primers respectively

    Returns two dictionaries:
        std_primer_dict_regex
            key: region+primer name
            val: primer sequence from 5' to 3'
        std_primer_dict
            key: region+primer name
            val: primer sequence from 5' to 3' for forward primers, 3' to 5' for reverse
    """

    std_primer_dict_regex = defaultdict(defaultdict)
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
                if line[0] == '>':
                    if 'R' in line: # If a primer is a reverse primer
                        rev_flag = True
                    key = line[1:]
                else:
                    if rev_flag:
                        rev_conv = str(Seq(line).reverse_complement())
                        line = rev_conv
                        rev_flag = False

                    primer = primer_regex_query_builder(line)
                    std_primer_dict_regex[region][key] = primer
                    std_primer_dict[region][key] = line

    return std_primer_dict_regex, std_primer_dict

def run_primer_agrep_once(input_path, input_primer, strand, mismatches=1):
    """
    Run the primer_agrep script.

    Takes one primer, strand, and fastq input
    Returns number of agrep matches
    """

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

def get_primer_props(std_primer_dict_regex, input_path):
    """
    Look for the standard primers in the input fastq file.

    Will loop through the dictionary of primers, using approximate grep (agrep) to find matching primers.
    If a std primer is present above a set threshold proportion, it is collected. Both strands are searched for.
    If there is an std primer for both the F and R strands, the maximum prop for each strand is chosen and the pair
    is output as a combination.

    Returns a list containing two elements:
        max_region: the amplified region the chosen primers belong to
        max_primers: dictionary containing the F and/or R primers that were chosen
    """

    threshold = 0.60 # Arbitrary threshold for collecting a matched primer
    read_count = get_read_count(input_path, 'fastq') # Get read count of fastq file to calculate proportion with
    res_dict = defaultdict(defaultdict)

    # Loop through every primer region
    for region, primer in std_primer_dict_regex.items():
        res_dict[region]['F'] = {}
        res_dict[region]['R'] = {}

        # Loop through every primer of a certain region
        for primer_name, primer_seq in primer.items():
            
            region_name_str = f'{region};{primer_name}'
            primer_count = 0.0

            if 'F' in primer_name:
                primer_count = run_primer_agrep_once(input_path, primer_seq, 'F', 1) # Get proportion of a F primer with agrep
            elif 'R' in primer_name:
                primer_count = run_primer_agrep_once(input_path, primer_seq, 'R', 1) # Get proportion of a R primer with agrep

            try:
                primer_prop = primer_count / read_count
            except ZeroDivisionError:
                primer_prop = 0

            if 'F' in primer_name:
                if primer_prop > threshold: # Only collect primer if it's above threshold
                    res_dict[region]['F'][primer_name] = primer_prop
            elif 'R' in primer_name:
                if primer_prop > threshold: # Only collect primer if it's above threshold
                    res_dict[region]['R'][primer_name] = primer_prop

            print(f'{region_name_str}: {primer_prop}')
        
        # If an F or/and R primer wasn't found then just remove it from the dictionary
        if res_dict[region]['F'] == {}:
            res_dict[region].pop('F')
        if res_dict[region]['R'] == {}:
            res_dict[region].pop('R')
        

    singles = defaultdict(str)
    doubles = defaultdict(list)

    double_status = False # Flag for whether primers were found on both strands

    #  Loop through every collected primer and put primers in singles or doubles
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
    
    # if at least one pair of primers was collected
    if double_status:
        for region in doubles: # Loop through all pairs of primers and choose the best one
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
        for region in singles: # Choose the best single primer
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



def save_out(results, sample_id, output, std_primer_dict):
    """
    Save found std primers into a fasta file.
    """

    with open(f'{output}/{sample_id}_std_primer_out.txt', 'w') as fw_out, open(f'{output}/{sample_id}_std_primers.fasta', 'w') as fw_seq:
        if results == []:
            fw_out.write(f'')
            fw_seq.write(f'')
        
        elif len(results) == 2:
            region = results[0]
            primer_name = list(results[1].keys())[0]
            primer_prop = results[1][list(results[1].keys())[0]]
            seq = std_primer_dict[region][primer_name]
            if 'R' in primer_name:
                seq = str(Seq(seq).reverse_complement())
            fw_out.write(f'{region}\n')
            fw_out.write(f'{primer_name}: {primer_prop}')

            fw_seq.write(f'>{primer_name}\n{seq}')
            
        
        elif len(results) == 3:
            region = results[0]
            f_primer_name = list(results[1].keys())[0]
            f_primer_prop = results[1][list(results[1].keys())[0]]
            f_seq = std_primer_dict[region][f_primer_name]
            r_primer_name = list(results[2].keys())[0]
            r_primer_prop = results[2][list(results[2].keys())[0]]
            r_seq = std_primer_dict[region][r_primer_name]
            r_seq = str(Seq(r_seq).reverse_complement())
            

            fw_out.write(f'{region}\n')
            fw_out.write(f'{f_primer_name}: {f_primer_prop}\n')
            fw_out.write(f'{r_primer_name}: {r_primer_prop}')

            fw_seq.write(f'>{f_primer_name}\n{f_seq}\n')
            fw_seq.write(f'>{r_primer_name}\n{r_seq}\n')

    
def main():
    
    _INPUT, _SAMPLE, _OUTPUT = parse_args()
    std_primer_dict_regex, std_primer_dict = parse_std_primers() # Parse std primer library into dictionaries
    results = get_primer_props(std_primer_dict_regex, _INPUT) # Find all the std primers in the input and select most common
    save_out(results, _SAMPLE, _OUTPUT, std_primer_dict)
    

if __name__ == "__main__":
    main()