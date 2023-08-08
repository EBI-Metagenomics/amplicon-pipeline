
include { fastp } from '../modules/fastp.nf'
include { seqprep_merge } from '../modules/seqprep_merge.nf'
include { fastq_to_fasta } from '../modules/fastq_to_fasta.nf'

workflow QC {

    take:
        project
        reads
        outdir

    main:
        fastp_input = project.combine(reads)
        fastp(fastp_input, outdir)
        seqprep_merge(fastp.out.cleaned_fastq, outdir)
        fastq_to_fasta(seqprep_merge.out.merged_fastq, outdir)
    
    emit:
        fastp_cleaned_fastq = fastp.out.cleaned_fastq
        merged_reads = seqprep_merge.out.merged_fastq
        merged_fasta = fastq_to_fasta.out.merged_fasta
    
}