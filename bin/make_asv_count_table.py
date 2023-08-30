
import argparse
from collections import defaultdict
import pandas as pd

_TAXA = "/hps/nobackup/rdf/metagenomics/service-team/users/chrisata/asv_nf_testing/merged/ERP122862/ERR4334396_taxa.tsv"
_FWD = "/hps/nobackup/rdf/metagenomics/service-team/users/chrisata/asv_nf_testing/merged/ERP122862/ERR4334396_1_map.txt"
_REV = "/hps/nobackup/rdf/metagenomics/service-team/users/chrisata/asv_nf_testing/merged/ERP122862/ERR4334396_2_map.txt"
_AMP = "/hps/nobackup/rdf/metagenomics/service-team/users/chrisata/nf-work/e7/caeed264657ef2ebdc6e9a86fa1ff5/ERR4334396_MERGED.V3-V4.txt"
_HEADERS = "/hps/software/users/rdf/metagenomics/service-team/users/chrisata/asv_gen/ERR4334396_1_headers.txt"
_SAMPLE = "ERR4334396"

# def parse_args():

#     parser = argparse.ArgumentParser()

#     parser.add_argument("-t", "--taxa", required=True, type=str, help="Path to merged FASTA to look for primers")
#     parser.add_argument("-s", "--sample", required=True, type=str, help="Sample ID")
#     parser.add_argument("-o", "--output", required=True, type=str, help="Output path")
#     args = parser.parse_args()
  
#     _INPUT = args.input
#     _SAMPLE = args.sample
#     _OUTPUT = args.output

#     return _INPUT, _SAMPLE, _OUTPUT


def main():
    
    taxa_df = pd.read_csv(_TAXA, sep="\t", dtype=str)
    taxa_df = taxa_df.fillna("0")
    taxa_df = taxa_df.sort_values(["Kingdom", "Phylum", "Class", "Order", "Genus"], ascending=True)

    amp_reads = [ read.strip()for read in list(open(_AMP, "r")) ]
    headers = [ read.split(" ")[0][1:] for read in list(open(_HEADERS, "r")) ]
    amp_region = _AMP.split(".")[1]

    asv_dict = defaultdict(int)

    with open(_FWD, "r") as fwd_fr, open(_REV, "r") as rev_fr:

        counter = -1
        for line_fwd, line_rev in zip(fwd_fr, rev_fr):

            counter += 1
            line_fwd = line_fwd.strip()
            line_rev = line_rev.strip()

            fwd_asvs = line_fwd.split(",")
            rev_asvs = line_rev.split(",")

            asv_intersection = list(set(fwd_asvs).intersection(rev_asvs))
            
            if len(asv_intersection) == 0:
                continue
            
            if len(asv_intersection) == 1 and asv_intersection[0] == "0":
                continue
            
            if headers[counter] in amp_reads:
                asv_dict[int(asv_intersection[0]) - 1] += 1
    
    tax_assignment_dict = defaultdict(int)

    for i in range(len(taxa_df)):
        
        sorted_index = taxa_df.index[i]
        asv_count = asv_dict[sorted_index]

        if asv_count == 0:
            continue

        k = taxa_df.loc[sorted_index, "Kingdom"]
        p = taxa_df.loc[sorted_index, "Phylum"]
        c = taxa_df.loc[sorted_index, "Class"]
        o = taxa_df.loc[sorted_index, "Order"]
        f = taxa_df.loc[sorted_index, "Family"]
        g = taxa_df.loc[sorted_index, "Genus"]        

        tax_assignment = ""

        while True:

            if k != "0":
                k = "_".join(k.split(" "))
                if k != "Archaea" and k != "Bacteria":
                    tax_assignment += f"sk__Eukaryota\tk__{k}"
                else:
                    tax_assignment += f"sk__{k}\tk__"
            else:
                break

            if p != "0":
                p = "_".join(p.split(" "))
                tax_assignment += f"\tp__{p}"
            else:
                break
            if c != "0":
                c = "_".join(c.split(" "))
                tax_assignment += f"\tc__{c}"
            else:
                break
            if o != "0":
                o = "_".join(o.split(" "))
                tax_assignment += f"\to__{o}"
            else:
                break
            if f != "0":
                f = "_".join(f.split(" "))
                tax_assignment += f"\tf__{f}"
            else:
                break
            if g != "0":
                g = "_".join(g.split(" "))
                tax_assignment += f"\tg__{g}"
            break

        tax_assignment_dict[tax_assignment] += asv_count



    # for asv_ind, asv_count in asv_dict.items():

    #     k = taxa_df.iloc[asv_ind, 1]
    #     p = taxa_df.iloc[asv_ind, 2]
    #     c = taxa_df.iloc[asv_ind, 3]
    #     o = taxa_df.iloc[asv_ind, 4]
    #     f = taxa_df.iloc[asv_ind, 5]
    #     g = taxa_df.iloc[asv_ind, 6]        

    #     tax_assignment = ""

    #     while True:

    #         if not pd.isna(k):
    #             k = "_".join(k.split(" "))
    #             if k != "Archaea" and k != "Bacteria":
    #                 tax_assignment += f"sk__Eukaryota\tk__{k}"
    #             else:
    #                 tax_assignment += f"sk__{k}\tk__"
    #         else:
    #             break

    #         if not pd.isna(p):
    #             p = "_".join(p.split(" "))
    #             tax_assignment += f"\tp__{p}"
    #         else:
    #             break
    #         if not pd.isna(c):
    #             c = "_".join(c.split(" "))
    #             tax_assignment += f"\tc__{c}"
    #         else:
    #             break
    #         if not pd.isna(o):
    #             o = "_".join(o.split(" "))
    #             tax_assignment += f"\to__{o}"
    #         else:
    #             break
    #         if not pd.isna(f):
    #             f = "_".join(f.split(" "))
    #             tax_assignment += f"\tf__{f}"
    #         else:
    #             break
    #         if not pd.isna(g):
    #             g = "_".join(g.split(" "))
    #             tax_assignment += f"\tg__{g}"
    #         break

    #     tax_assignment_dict[tax_assignment] += asv_count

    with open(f"./{_SAMPLE}_{amp_region}_asv_krona_counts.txt", "w") as fw:
        for tax_assignment, count in tax_assignment_dict.items():
            fw.write(f"{count}\t{tax_assignment}\n")
    
        # print(taxa_df.iloc[asv_ind, 0])


if __name__ == "__main__":
    main()