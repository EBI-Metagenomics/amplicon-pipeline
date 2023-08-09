nextflow.enable.dsl=2

// nextflow run nextflow/primer_trimming.nf --path /hps/nobackup/rdf/metagenomics/service-team/users/chrisata/asv_datasets/ERP123542/raw --project ERP123542 --outdir /hps/nobackup/rdf/metagenomics/service-team/users/chrisata/asv_nf_testing/merged
// nextflow run nextflow/primer_trimming.nf -profile lsf --path /hps/nobackup/rdf/metagenomics/service-team/users/chrisata/asv_datasets/ERP123542/raw --project ERP123542 --outdir /hps/nobackup/rdf/metagenomics/service-team/users/chrisata/asv_nf_testing

include { cutadapt } from './modules/cutadapt.nf'
include { concat_primers } from './modules/concat_primers.nf'

include { QC } from './subworkflows/qc_swf.nf'
include { PRIMER_IDENTIFICATION } from './subworkflows/primer_identification_swf.nf'
include { AUTOMATIC_PRIMER_PREDICTION } from './subworkflows/automatic_primer_trimming.nf'
include { CMSEARCH } from './subworkflows/cmsearch_swf.nf'

params.path = null
params.project = null
params.outdir = null

reads = Channel.fromFilePairs( "${params.path}/ERR447178*_{1,2}.fastq.gz" )
project = Channel.value( params.project )
outdir = params.outdir

workflow {

    // Workflow for running the primer-trimming survey on sequencing reads

    // Quality control
    // TODO: need to get this working for single-end reads
    QC(
        project,
        reads,
        outdir
    )

    // cmsearch to get rRNA matches
    CMSEARCH(
        QC.out.merged_fasta,
        outdir
    )

    // CMSEARCH.out.cmsearch_deoverlap_out.view()

    // TODO: amplified region process will go here

    // Automatic primer identification
    PRIMER_IDENTIFICATION(
        QC.out.merged_reads,
        outdir
    )

    // Join flags with merged fastq files
    auto_trimming_input = PRIMER_IDENTIFICATION.out.conductor_out
    .join(QC.out.merged_reads)


    // Run subworkflow for automatic primer trimming
    // Outputs either empty fasta file or fasta file containing predicted primers
    AUTOMATIC_PRIMER_TRIMMING(
        auto_trimming_input,
        outdir
        )

    // Join auto flags to std flags and generated a concatenated fasta file containing primers to trim off
    // This can contain any valid combination of stranded std/auto primers 
    concat_input = PRIMER_IDENTIFICATION.out.std_primer_out
    .join(AUTOMATIC_PRIMER_TRIMMING.out.auto_primer_trimming_out)
    concat_primers(concat_input, outdir)

    // Join concatenated primers to the fastp-cleaned paired reads files and run cutadapt on them
    cutadapt_input = concat_primers.out.concat_primers_out
    .join(QC.out.fastp_cleaned_fastq)
    cutadapt(cutadapt_input, outdir)

    // Just some logging for myself, will delete this eventually
    // final_out = fastp.out.cleaned_fastq
    // .join(cutadapt.out.cutadapt_out, remainder: true)
    // .map( { if (it[4] == null) { tuple(it[0], it[1], it[2], it[3]) } else { tuple(it[0], it[1], it[4], it[5]) }} )


}