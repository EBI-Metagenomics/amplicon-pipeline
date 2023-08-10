
include { fastp } from '../modules/fastp.nf'
include { seqprep_merge } from '../modules/seqprep_merge.nf'
include { trimmomatic_SE } from '../modules/trimmomatic_SE.nf'
include { fastq_to_fasta } from '../modules/fastq_to_fasta.nf'

workflow QC {

    // Quality control subworkflow
    // Removes adapters with fastp, merges paired-end reads with seqprep, converts resulting fastq file to fasta

    take:
        project
        reads
        outdir

    main:
        fastp_input = project.combine(reads)
        fastp(fastp_input, outdir)
        seqprep_merge(fastp.out.cleaned_fastq, outdir)
        trimmomatic_SE(seqprep_merge.out.merged_fastq, outdir)
        fastq_to_fasta(trimmomatic_SE.out.trimmed_fastq, outdir)
    
    emit:
        fastp_cleaned_fastq = fastp.out.cleaned_fastq
        merged_reads = trimmomatic_SE.out.trimmed_fastq
        merged_fasta = fastq_to_fasta.out.merged_fasta
    
}