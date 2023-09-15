
include { FASTP } from '../../modules/nf-core/fastp/main.nf'
include { SEQPREP_MERGE } from '../../modules/local/seqprep_merge.nf'
include { TRIMMOMATIC_SE } from '../../modules/local/trimmomatic_SE.nf'
include { FASTQ_TO_FASTA } from '../../modules/local/fastq_to_fasta.nf'

workflow QC {

    // Quality control subworkflow
    // Removes adapters with fastp, merges paired-end reads with seqprep, converts resulting fastq file to fasta

    take:
        project
        reads
        outdir

    main:
        // fastp_input = project
        //               .combine(reads)

        // fastp_input.view()

        FASTP(
            reads,
        )

        // SEQPREP_MERGE(
        //     FASTP.out.cleaned_fastq,
        //     outdir
        // )

        // TRIMMOMATIC_SE(
        //     SEQPREP_MERGE.out.merged_fastq, 
        //     outdir
        // )

        // FASTQ_TO_FASTA(
        //     TRIMMOMATIC_SE.out.trimmed_fastq,
        //     outdir
        // )
    
    // emit:
        // fastp_cleaned_fastq = FASTP.out.cleaned_fastq
        // merged_reads = TRIMMOMATIC_SE.out.trimmed_fastq
        // merged_fasta = FASTQ_TO_FASTA.out.merged_fasta
    
}