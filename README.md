# asv-gen

Extremely-in-progress pipeline for the generation of Amplicon Sequence Variants from amplicon data. Will write a far more detailed README in the near future when it's closer to a good first draft.

At the moment the pipeline mainly consists of the primer workflow, which is mostly complete and functioning. See below for instructions for running it.

## Requirements

At the moment the prerequisites are Singularity and a micromamba environment located at:

`/hps/software/users/rdf/metagenomics/service-team/software/micromamba/envs/asv-test`

## How to run

At the moment you can run the primer trimming pipeline like this:

`nextflow run nextflow/primer_trimming.nf --path {directory/containing/fastq} --project {project_accession} --outdir {path/to/out}`

