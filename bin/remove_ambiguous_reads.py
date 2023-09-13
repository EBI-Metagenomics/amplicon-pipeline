
# Script removes any reads with ambiguous bases (Ns) for the purpose of DADA2

import argparse
import gzip

from Bio import SeqIO, bgzf

_INPUT1 = "/hps/nobackup/rdf/metagenomics/service-team/users/chrisata/nf-work/f8/570e8803fa29a407fa1786206f5ca2/SRR17062740_1.cutadapt.fastq.gz"
_INPUT2 = "/hps/nobackup/rdf/metagenomics/service-team/users/chrisata/nf-work/f8/570e8803fa29a407fa1786206f5ca2/SRR17062740_2.cutadapt.fastq.gz"


def parse_args():

    parser = argparse.ArgumentParser()

    parser.add_argument("-f", "--fwd", required=True, type=str, help="Path to fastq file to check for primers")
    parser.add_argument("-r", "--rev", required=True, type=str, help="Fasta")
    parser.add_argument("-s", "--sample", required=True, type=str, help="Sample ID")
    args = parser.parse_args()
    
    _FWD = args.fwd
    _REV = args.rev
    _SAMPLE = args.sample

    return _FWD, _REV, _SAMPLE


def main():

    _FWD, _REV, _SAMPLE = parse_args()
    
    remove_lst = []

    with gzip.open(_FWD, "rt") as fwd_handle, gzip.open(_REV, "rt") as rev_handle:
        fwd_reads = SeqIO.to_dict(SeqIO.parse(fwd_handle, "fastq"))
        rev_reads = SeqIO.to_dict(SeqIO.parse(rev_handle, "fastq"))

        for read_id in fwd_reads.keys():

            if "N" in str(fwd_reads[read_id].seq):
                print(read_id)
                remove_lst.append(read_id)
                continue
            elif "N" in str(rev_reads[read_id].seq):
                print(read_id)
                remove_lst.append(read_id)
                continue


    [ fwd_reads.pop(read_id) for read_id in remove_lst ]
    [ rev_reads.pop(read_id) for read_id in remove_lst ]
    
    with bgzf.BgzfWriter(f"./{_SAMPLE}_noambig_1.fastq.gz", "wb") as fwd_handle, bgzf.BgzfWriter(f"./{_SAMPLE}_noambig_2.fastq.gz", "wb") as rev_handle:
        SeqIO.write(sequences=fwd_reads.values(), handle=fwd_handle, format="fastq")
        SeqIO.write(sequences=rev_reads.values(), handle=rev_handle, format="fastq")

if __name__ == "__main__":
    main()