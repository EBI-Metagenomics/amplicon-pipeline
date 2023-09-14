# asv-gen

This repository contains the in-progress Nextflow rework of the V5 MGnify amplicon annotation pipeline, which will include the addition of Amplicon Sequence Variant (ASV) annotation as a new major feature.

The reworked pipeline can be divided into seven different subworkflows, four of which already exist in the current V5 amplicon pipeline:

* Quality Control (V5)
* rRNA Prediction and SSU+LSU extraction (V5)
* ITS extraction (V5)
* SSU+LSU+ITS taxonomic classification and visualisation (V5)

The new subworkflows consist of:

* Amplified region inference
* Primer trimming
* ASV generation, classification, and visualisation

![V6 Schema](https://github.com/EBI-Metagenomics/asv-gen/assets/34323164/4b42b6c0-dfd2-4fd7-a04f-942a1f3904fa)

## Completion

The different subworkflows' levels of completion:

- Quality Control (V5) :white_check_mark:
- rRNA Prediction and SSU+LSU extraction (V5) :white_check_mark:
- ITS extraction (V5) :white_check_mark:
- SSU+LSU+ITS taxonomic classification and visualisation (V5) :white_check_mark:
- Amplified region inference :white_check_mark:
- Primer trimming :white_check_mark:
- ASV generation, classification, and visualisation :white_check_mark:


## Requirements

At the moment the prerequisites are Singularity and a micromamba environment located at:

`/hps/software/users/rdf/metagenomics/service-team/software/micromamba/envs/asv-test`


## How to run

You can run the current version of the pipeline like this:

`nextflow run main.nf --path {directory/containing/fastq} --project {project_accession} --outdir {path/to/out}`

## Output Directory Structure

Example output directory structure for one run:

```
.
├── amplified-region-inference
│   ├── ERR4334351.trimmed.fastq_16S-V3-V4_extracted.fastq.gz
│   └── ERR4334351.tsv
├── asv-gen
│   ├── 16S-V3-V4
│   │   └── ERR4334351_16S.V3-V4_asv_krona_counts.txt
│   ├── ERR4334351_proportion_chimeric.txt
│   ├── ERR4334351_proportion_matched.txt
│   └── ERR4334351_taxa.tsv
├── primer-identification
│   ├── ERR4334351_16S-V3-V4_auto_primers.fasta
│   ├── ERR4334351_16S-V3-V4_std_primer_out.txt
│   ├── ERR4334351_16S-V3-V4_std_primers.fasta
│   ├── ERR4334351_1.cutadapt.fastq.gz
│   ├── ERR4334351_2.cutadapt.fastq.gz
│   ├── ERR4334351_final_concat_primers.fasta
│   ├── ERR4334351_primer_validation.tsv
│   └── ERR4334351_trimming_conductor_out_16S-V3-V4.txt
├── QC
│   ├── ERR4334351.fasta
│   └── ERR4334351.trimmed.fastq.gz
├── sequence-categorisation
│   ├── ERR4334351.cmsearch_matches.tbl.deoverlapped
│   └── ERR4334351_SSU.fasta
└── taxonomy-summary
    ├── DADA2-SILVA
    │   └── ERR4334351_16S.V3-V4_asv_krona_counts_krona.html
    ├── ITSonedb
    │   ├── ERR4334351_ITS_masked_krona.html
    │   ├── ERR4334351_ITS_masked.notaxid.tsv
    │   ├── ERR4334351_ITS_masked.tsv
    │   └── ERR4334351_ITS_masked.txt
    ├── SSU
    │   ├── ERR4334351_SSU_krona.html
    │   ├── ERR4334351_SSU.notaxid.tsv
    │   ├── ERR4334351_SSU.tsv
    │   └── ERR4334351_SSU.txt
    └── UNITE
        ├── ERR4334351_ITS_masked_krona.html
        ├── ERR4334351_ITS_masked.notaxid.tsv
        ├── ERR4334351_ITS_masked.tsv
        └── ERR4334351_ITS_masked.txt
```
