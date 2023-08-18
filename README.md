# asv-gen

This repository contains the extremely-in-progress Nextflow rework of the V5 MGnify amplicon annotation pipeline, which will include the addition of Amplicon Sequence Variant annotation as a new major feature.

The reworked pipeline can be divided into seven different subworkflows, four of which already exist in the current V5 amplicon pipeline:

* Quality Control (V5)
* rRNA Prediction and SSU+LSU extraction (V5)
* ITS extraction (V5)
* SSU+LSU+ITS taxonomic classification and visualisation (V5)

The new subworkflows will consist of:

* Amplified region inference
* Primer trimming
* ASV generation, classification, and visualisation

![V6 Schema](https://github.com/EBI-Metagenomics/asv-gen/assets/34323164/a6ef22eb-7967-4d1c-b635-0468eb11e174)

## Completion

The different subworkflow's levels of completion 

- Quality Control (V5) :white_check_mark:
- rRNA Prediction and SSU+LSU extraction (V5) :white_check_mark:
- ITS extraction (V5) :x:
- SSU+LSU+ITS taxonomic classification and visualisation (V5) :white_check_mark:
- Amplified region inference ✅
- Primer trimming :soon:
- ASV generation, classification, and visualisation :soon:

## Requirements

At the moment the prerequisites are Singularity and a micromamba environment located at:

`/hps/software/users/rdf/metagenomics/service-team/software/micromamba/envs/asv-test`

## How to run

At the moment you can run the current version of the pipeline like this:

`nextflow run nextflow/primer_trimming.nf --path {directory/containing/fastq} --project {project_accession} --outdir {path/to/out}`

