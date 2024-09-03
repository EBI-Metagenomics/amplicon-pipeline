# MGnify amplicon analysis pipeline

This repository contains the **beta version** of the v6.0 [MGnify](https://www.ebi.ac.uk/metagenomics) amplicon analysis pipeline. It is, first and foremost, a refactor of the existing [v5.0 amplicon analysis pipeline](https://github.com/EBI-Metagenomics/pipeline-v5), replacing CWL with [Nextflow](https://www.nextflow.io/) as its workflow management system. This pipeline re-implements all [existing closed-reference v5.0 features](https://docs.mgnify.org/src/docs/analysis.html#amplicon-analysis-pipeline), and makes multiple significant changes and additions.

![V6 Schema](assets/v6_amplicon_schema.png)

## Pipeline description

### Features

The amplicon analysis pipeline v6.0 re-implements all of the existing features from v5.0:

- Reads quality control
- rRNA sequence extraction using [Infernal/cmsearch](https://github.com/EddyRivasLab/infernal/tree/master)
- Closed-reference-based taxonomic classification and visualiation of rRNA using [MAPseq](https://github.com/meringlab/MAPseq) and [Krona](https://github.com/marbl/Krona)

The amplicon analysis pipeline v6.0 also contains multiple significant changes:

- Refactoring from CWL to [Nextflow](https://www.nextflow.io/) for pipeline definition
- Simplification the reads quality control using [fastp](https://github.com/OpenGene/fastp)
- Automatic amplified region inference for 16S and 18S rRNA
- Automatic primer identification, trimming, and validation
- Addition of Amplicon Sequence Variant (ASV) calling using [DADA2](https://benjjneb.github.io/dada2/index.html)
- Taxonomic classification and visualisation of ASVs using [MAPseq](https://github.com/meringlab/MAPseq) and [Krona](https://github.com/marbl/Krona) to complement the existing closed-reference analysis
- Addition of [PR2](https://pr2-database.org/) as a reference database
- Updating of existing reference databases ([SILVA](https://www.arb-silva.de/), [UNITE](https://unite.ut.ee/), [ITSoneDB](https://itsonedb.cloud.ba.infn.it), [Rfam](https://rfam.org/))

### Valid amplicons

At this stage, the only sequence amplicons that this pipeline is built for are:

| Amplicon 	| Closed-reference analysis 	| ASV analysis 	|
|:--------:	|:-------------------------:	|:------------:	|
|    16S   	|             ✓             	|       ✓      	|
|    18S   	|             ✓             	|       ✓      	|
|    LSU   	|             ✓             	|       ✗      	|
|    ITS   	|             ✓             	|       ✗      	|


### Tools

| Tool                      	| Version  	| Purpose                                                	|
|---------------------------	|----------	|--------------------------------------------------------	|
| [fastp](https://github.com/OpenGene/fastp)                     	| 0.23.4   	| Read quality control                                   	|
| [seqtk](https://github.com/lh3/seqtk)                     	| 1.3-r106 	| FASTQ file manipulation                                	|
| [easel](https://github.com/EddyRivasLab/easel)                     	| 0.49     	| FASTA file manipulation                                	|
| [bedtools](https://bedtools.readthedocs.io/en/latest/)                  	| 2.30.0   	| FASTA sequence masking                                 	|
| [Infernal/cmsearch](https://github.com/EddyRivasLab/infernal/tree/master)                  	| 1.1.5    	| rRNA sequence searching                                	|
| [cmsearch_tblout_deoverlap](https://github.com/nawrockie/cmsearch_tblout_deoverlap/tree/master) 	| 0.09     	| Deoverlapping of cmsearch results                      	|
| [MAPseq](https://github.com/meringlab/MAPseq)                    	| 2.1.1b   	| Reference-based taxonomic classification of rRNA       	|
| [Krona](https://github.com/marbl/Krona)                    	| 2.8.1    	| Krona chart visualisation                              	|
| [cutadapt](https://cutadapt.readthedocs.io/en/stable/)                  	| 4.6      	| Primer trimming                                        	|
| [R](https://www.r-project.org/)                         	| 4.3.3    	| R programming language (runs DADA2)                    	|
| [DADA2](https://benjjneb.github.io/dada2/index.html)                     	| 1.30.0   	| ASV calling                                            	|
| [mgnify-pipelines-toolkit](https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit)  	| 0.1.5    	| Toolkit containing various in-house processing scripts 	|

### Reference databases

This pipeline uses five different reference databases. The files the pipeline uses are processed from the raw files available on each database's website, for use by MAPseq and cmsearch. We provide ready-made versions of these processed files on our FTP, which you can find here:

| Reference database 	| Version 	| Purpose                               	| Processed file paths                                                                                                                                          	|
|--------------------	|---------	|---------------------------------------	|---------------------------------------------------------------------------------------------------------------------------------------------------------------	|
| [SILVA](https://www.arb-silva.de/)              	| 138.1   	| 16S+18S+LSU rRNA database                 	| https://ftp.ebi.ac.uk/pub/databases/metagenomics/pipelines/tool-dbs/silva-ssu/ https://ftp.ebi.ac.uk/pub/databases/metagenomics/pipelines/tool-dbs/silva-lsu/ 	|
| [PR2](https://pr2-database.org/)                	| 5.0     	| Protist-focused 18S+16S rRNA database 	| https://ftp.ebi.ac.uk/pub/databases/metagenomics/pipelines/tool-dbs/pr2/                                                                                      	|
| [UNITE](https://unite.ut.ee/)              	| 9.0     	| ITS database                          	| https://ftp.ebi.ac.uk/pub/databases/metagenomics/pipelines/tool-dbs/unite/                                                                                    	|
| [ITSoneDB](https://itsonedb.cloud.ba.infn.it)           	| 1.141   	| ITS database                          	| https://ftp.ebi.ac.uk/pub/databases/metagenomics/pipelines/tool-dbs/itsonedb/                                                                                 	|
| [Rfam](https://rfam.org/)               	| 14.10   	| RNA family profile database           	| https://ftp.ebi.ac.uk/pub/databases/metagenomics/pipelines/tool-dbs/rfam/                                                                                     	|

## How to run

### Requirements

At the moment the only prerequisites for running it are Nextflow and [Docker](https://www.docker.com/)/[Singularity](https://docs.sylabs.io/guides/3.5/user-guide/introduction.html), since all of the Nextflow processes use pre-built containers.

### Input shape

The input data for the pipeline is amplicon sequencing reads (either paired-end or single-end) in the form of FASTQ files. These files should be specified using a `.csv` samplesheet file with this format:

```
sample,fastq_1,fastq_2,single_end
SRR9674618,/path/to/reads/SRR9674618.fastq.gz,,true
SRR17062740,/path/to/reads/SRR17062740_1.fastq.gz,/path/to/reads/SRR17062740_2.fastq.gz,false
```

### Execution

You can run the current version of the pipeline on SLURM like this:

`nextflow run -profile codon_slurm main.nf --input /path/to/samplesheet.csv --outdir /path/to/outdir main.nf`

### Output directory structure

Example output directory structure for one run:

```
├── amplified-region-inference
│   ├── ERR2356289_16S-V4-V5_extracted.fastq.gz
│   ├── ERR2356289.16S.V4-V5.txt
│   └── ERR2356289.tsv
├── asv
│   ├── 16S-V4-V5
│   │   ├── ERR2356289_16S.V4-V5_asv_read_counts.tsv
│   │   ├── ERR2356289_16S.V4-V5_pr2_asv_krona_counts.txt
│   │   └── ERR2356289_16S.V4-V5_silva_asv_krona_counts.txt
│   ├── ERR2356289_asv_seqs.fasta
│   ├── ERR2356289_DADA2-PR2_asv_taxa.tsv
│   ├── ERR2356289_DADA2-SILVA_asv_taxa.tsv
│   ├── ERR2356289_proportion_chimeric.txt
│   └── ERR2356289_proportion_matched.txt
├── primer-identification
│   ├── ERR2356289_16S-V4-V5_auto_primers.fasta
│   ├── ERR2356289_16S-V4-V5_std_primers.fasta
│   ├── ERR2356289_1.trim.fastq.gz
│   ├── ERR2356289_2.trim.fastq.gz
│   ├── ERR2356289_final_concat_primers.fasta
│   ├── ERR2356289_rev_comp_se_primers.fasta
│   └── ERR2356289_trimming_conductor_out_16S-V4-V5.txt
├── qc
│   ├── ERR2356289.merged.fastq.gz
│   └── ERR2356289.seqtk-seq.fasta.gz
├── sequence-categorisation
│   ├── ERR2356289_SSU.fasta
│   └── ERR2356289.tblout.deoverlapped
└── taxonomy-summary
    ├── DADA2-PR2
    │   └── ERR2356289_16S-V4-V5.html
    ├── DADA2-SILVA
    │   └── ERR2356289_16S-V4-V5.html
    ├── ITSoneDB
    │   ├── ERR2356289.html
    │   └── ERR2356289.txt
    ├── PR2
    │   ├── ERR2356289.html
    │   └── ERR2356289.txt
    ├── SSU
    │   ├── ERR2356289.html
    │   └── ERR2356289.txt
    └── UNITE
        ├── ERR2356289.html
        └── ERR2356289.txt
```
