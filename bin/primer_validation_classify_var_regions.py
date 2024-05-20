#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright 2024 EMBL - European Bioinformatics Institute
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import argparse
from collections import defaultdict
import re

from Bio import SeqIO
import pandas as pd


regions_16S_bacteria = {
    "V1": [69, 92],
    "V2": [131, 239],
    "V3": [430, 487],
    "V4": [566, 672],
    "V5": [812, 869],
    "V6": [976, 1033],
    "V7": [1107, 1164],
    "V8": [1234, 1285],
    "V9": [1426, 1456],
}

regions_16S_archaea = {
    "V1": [61, 79],
    "V2": [114, 223],
    "V3": [397, 436],
    "V4": [516, 623],
    "V5": [763, 824],
    "V6": [932, 982],
    "V7": [1056, 1119],
    "V8": [1189, 1240],
    "V9": [1372, 1410],
}

regions_18S = {
    "V1": [69, 109],
    "V2": [136, 298],
    "V3": [474, 545],
    "V4": [627, 873],
    "V5": [1059, 1102],
    "V7": [1366, 1454],
    "V8": [1526, 1608],
    "V9": [1728, 1795],
}


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--input", required=True, type=str, help="Path to cmsearch_deoverlap_tblout file")
    parser.add_argument("-f", "--fasta", required=True, type=str, help="Path to concatenated primers fasta file")
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

            if rfam == "RF00177":
                gene = "16S"
                model = regions_16S_bacteria
            elif rfam == "RF01959":
                gene = "16S"
                model = regions_16S_archaea
            elif rfam == "RF01960":
                gene = "18S"
                model = regions_18S
            else:
                continue

            res_dict["Run"].append(_SAMPLE)
            res_dict["AssertionEvidence"].append("ECO_0000363")
            res_dict["AssertionMethod"].append("automatic assertion")

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
