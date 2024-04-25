# asv-gen

This repository contains the in-progress Nextflow rework of the V5 MGnify amplicon annotation pipeline, which will include the addition of Amplicon Sequence Variant (ASV) annotation as a new major feature under the V6 analysis pipelines.

The reworked pipeline can be divided into seven different subworkflows, four of which already exist in the current V5 amplicon pipeline:

- Quality Control (V5)
- rRNA Prediction and SSU+LSU extraction (V5)
- ITS extraction (V5)
- SSU+LSU+ITS taxonomic classification and visualisation (V5)

The new subworkflows consist of:

- Amplified region inference
- Primer trimming
- ASV generation, classification, and visualisation

![V6 Schema](https://github.com/EBI-Metagenomics/asv-gen/assets/34323164/a6ef22eb-7967-4d1c-b635-0468eb11e174)

## Requirements

At the moment the prerequisites are Singularity and a micromamba environment located at:

`/hps/software/users/rdf/metagenomics/service-team/software/micromamba/envs/asv-test`

## Input shape

The input data to the pipeline should be a `samplesheet.csv` file with this format:

```
sample,fastq_1,fastq_2,single_end
SRR9674618,/path/to/reads/SRR9674618.fastq.gz,,true
SRR17062740,/path/to/reads/SRR17062740_1.fastq.gz,/path/to/reads/SRR17062740_2.fastq.gz,false
```

## How to run

You can run the current version of the pipeline on SLURM like this:

`nextflow run -profile codon_slurm main.nf --input /path/to/samplesheet --outdir /path/to/outdir`

## Output directory structure

Example output directory structure for one run:

```
.
├── amplified-region-inference
│   ├── ERR4334351_16S-V3-V4_extracted.fastq.gz
│   ├── ERR4334351.16S.V3-V4.txt
│   └── ERR4334351.tsv
├── asv-gen
│   ├── 16S-V3-V4
│   ├── ERR4334351_asvs.fasta
│   ├── ERR4334351_proportion_chimeric.txt
│   └── ERR4334351_proportion_matched.txt
├── primer-identification
│   ├── ERR4334351_16S-V3-V4_auto_primers.fasta
│   ├── ERR4334351_16S-V3-V4_std_primers.fasta
│   ├── ERR4334351_1.trim.fastq.gz
│   ├── ERR4334351_2.trim.fastq.gz
│   ├── ERR4334351_final_concat_primers.fasta
│   ├── ERR4334351_primer_validation.tsv
│   ├── ERR4334351_rev_comp_se_primers.fasta
│   └── ERR4334351_trimming_conductor_out_16S-V3-V4.txt
├── qc
│   ├── ERR4334351.merged.fastq.gz
│   └── ERR4334351.seqtk-seq.fasta.gz
├── sequence-categorisation
│   ├── ERR4334351_SSU.fasta
│   └── ERR4334351.tblout.deoverlapped
└── taxonomy-summary
    ├── DADA2-PR2
    ├── DADA2-SILVA
    ├── ITSonedb
    ├── PR2
    ├── SSU
    └── UNITE
```
