
import argparse
from collections import defaultdict

import pandas as pd

def parse_args():

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True, type=str, help="Input")
    parser.add_argument("-l", "--label", choices=['DADA2-SILVA', 'DADA2-PR2'], required=True, type=str, help="Database label - either DADA2-SILVA or DADA2-PR2")
    parser.add_argument("-s", "--sample", required=True, type=str, help="Sample ID")

    args = parser.parse_args()

    _INPUT = args.input
    _LABEL = args.label
    _SAMPLE = args.sample

    return _INPUT, _LABEL, _SAMPLE

def parse_label(label):

    silva_short_ranks = ["sk__", "k__", "p__", "c__", "o__", "f__", "g__", "s__"]
    pr2_short_ranks = ["d__", "sg__", "dv__", "sdv__", "c__", "o__", "f__", "g__", "s__"]

    silva_long_ranks = ["Superkingdom", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
    pr2_long_ranks = ["Domain", "Supergroup", "Division", "Subdivision", "Class", "Order", "Family", "Genus", "Species"]

    chosen_short_ranks = ''
    chosen_long_ranks = ''

    if label == 'DADA2-SILVA':
        chosen_short_ranks = silva_short_ranks
        chosen_long_ranks = silva_long_ranks
    elif label == 'DADA2-PR2':
        chosen_short_ranks = pr2_short_ranks
        chosen_long_ranks = pr2_long_ranks
    else:
        print("Incorrect database label - exiting")
        exit(1)

    return chosen_short_ranks, chosen_long_ranks

def parse_mapseq(mseq_df, short_ranks, long_ranks):

    res_dict = defaultdict(list)

    for i in range(len(mseq_df)):
        asv_id = mseq_df.iloc[i, 0]
        tax_ass = mseq_df.iloc[i, 1].split(';')

        res_dict['ASV'].append(asv_id)
        
        for j in range(len(short_ranks)):

            curr_rank = long_ranks[j]
            
            if j >= len(tax_ass):
                # This would only be true if the assigned taxonomy is shorter than the total reference database taxononmy
                # so fill each remaining rank with its respective short rank blank
                curr_tax = short_ranks[j]
            else:
                curr_tax = tax_ass[j]
            
            res_dict[curr_rank].append(curr_tax)
    res_df = pd.DataFrame.from_dict(res_dict)

    return(res_df)

def process_blank_tax_ends(res_df, ranks):
    # Necessary function as we want to replace consecutive blank assignments that start at the last rank as NAs
    # while avoiding making blanks in the middle as NAs

    for i in range(len(res_df)):
        last_empty_rank = ''
        currently_empty = False
        for j in reversed(range(len(ranks))): # Parse an assignment backwards, from Species all the way to Superkingdom/Domain
            curr_rank = res_df.iloc[i, j+1]
            if curr_rank in ranks:
                if last_empty_rank == '': # Last rank is empty, start window of consecutive blanks 
                    last_empty_rank = j+1
                    currently_empty = True
                elif currently_empty: # If we're in a window of consecutive blank assignments that started at the beginning
                    last_empty_rank = j+1
                else:
                    break
            else:
                break
        if last_empty_rank != '':
            res_df.iloc[i, last_empty_rank:] = 'NA'

    return res_df

def main():
    
    _INPUT, _LABEL, _SAMPLE = parse_args()

    mseq_df = pd.read_csv(_INPUT, header=1, sep='\t', usecols=[0, 13])

    short_ranks, long_ranks = parse_label(_LABEL)
    res_df = parse_mapseq(mseq_df, short_ranks, long_ranks)
    final_res_df = process_blank_tax_ends(res_df, short_ranks)

    final_res_df.to_csv(f"./{_SAMPLE}_{_LABEL}_asv_taxa.tsv", sep="\t", index=False)

if __name__ == "__main__":
    main()