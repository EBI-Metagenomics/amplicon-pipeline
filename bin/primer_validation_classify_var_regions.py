
import argparse
from collections import defaultdict
import re

from Bio import SeqIO
import pandas as pd


# _INPUT = "/hps/nobackup/rdf/metagenomics/service-team/users/chrisata/nf-work/36/6397218da76401d68fee87b15734cc/mock_V3-V4_V9_double_amp.cmsearch_matches.tbl.deoverlapped"
# _FASTA = "/hps/nobackup/rdf/metagenomics/service-team/users/chrisata/nf-work/5c/ce44ee574ce7acba820a5764c52d89/mock_V3-V4_V9_double_amp_final_concat_primers.fasta"
# _SAMPLE = "mock_V3-V4_V9_double_amp"

regions_16S_bacteria = {
    'V1': [69, 92],
    'V2': [131, 239],
    'V3': [430, 487],
    'V4': [566, 672],
    'V5': [812, 869],
    'V6': [976, 1033],
    'V7': [1107, 1164],
    'V8': [1234, 1285],
    'V9': [1426, 1456]
}

regions_16S_archaea = {
    'V1': [61, 79],
    'V2': [114, 223],
    'V3': [397, 436],
    'V4': [516, 623],
    'V5': [763, 824],
    'V6': [932, 982],
    'V7': [1056, 1119],
    'V8': [1189, 1240],
    'V9': [1372, 1410]
}

regions_18S = {
    'V1': [69, 109],
    'V2': [136, 298],
    'V3': [474, 545],
    'V4': [627, 873],
    'V5': [1059, 1102],
    'V7': [1366, 1454],
    'V8': [1526, 1608],
    'V9': [1728, 1795]
}


def parse_args():

    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--input", required=True, type=str, help="Path to fastq file to check for primers")
    parser.add_argument("-f", "--fasta", required=True, type=str, help="Fasta")
    parser.add_argument("-s", "--sample", required=True, type=str, help="Sample ID")
    args = parser.parse_args()
    
    _INPUT = args.input
    _FASTA = args.fasta
    _SAMPLE = args.sample

    return _INPUT, _FASTA, _SAMPLE


def get_amp_region(beg, end, strand, model):

    prev_region = ""

    for region, region_coords in model.items():

        region_beg = region_coords[0]
        region_end = region_coords[1]

        beg_diff = region_beg - beg
        end_diff = region_end - end

        if strand == "fwd":
            if beg_diff > 0:
                return region
        else:
            if beg_diff > 0:
                return prev_region
        
        prev_region = region

    return prev_region

def main():

    _INPUT, _FASTA, _SAMPLE = parse_args()    
    res_dict = defaultdict(list)
    fasta_dict = SeqIO.to_dict(SeqIO.parse(_FASTA, "fasta"))

    with open(_INPUT, "r") as fr:
        for line in fr:
            line = line.strip()            
            line = re.sub("[ \t]+", "\t", line)
            line_lst = line.split("\t")

            primer_name = line_lst[0]
            rfam = line_lst[3]
            beg = float(line_lst[5])
            end = float(line_lst[6])

            res_dict["Run"].append(_SAMPLE)
            res_dict["AssertionEvidence"].append("ECO_0000363")
            res_dict["AssertionMethod"].append("automatic assertion")

            if rfam == 'RF00177':
                gene = "16S"
                model = regions_16S_bacteria
            elif rfam == 'RF01959':
                gene = "16S"
                model = regions_16S_archaea
            elif rfam == 'RF01960':
                gene = "18S"
                model = regions_18S

            strand = ""

            if "F" in primer_name:
                strand = "fwd"
            elif "R" in primer_name: 
                strand = "rev"

            amp_region = get_amp_region(beg, end, strand, model)
            primer_seq = str(fasta_dict[primer_name].seq)

            res_dict["Gene"].append(gene)
            res_dict["VariableRegion"].append(amp_region)
            res_dict["PrimerName"].append(primer_name)
            res_dict["PrimerStrand"].append(strand)
            res_dict["PrimerSeq"].append(primer_seq)


    res_df = pd.DataFrame.from_dict(res_dict)
    res_df.to_csv(f"./{_SAMPLE}_primer_validation.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()