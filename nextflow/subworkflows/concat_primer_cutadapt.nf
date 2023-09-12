
include { CONCAT_PRIMERS } from '../modules/concat_primers.nf'
include { FINAL_CONCAT_PRIMERS } from '../modules/final_concat_primers.nf'
include { CUTADAPT } from '../modules/cutadapt.nf'

workflow CONCAT_PRIMER_CUTADAPT {
    
    take:
        concat_input
        fastp_cleaned_fastq
        outdir

    main:
        CONCAT_PRIMERS(
            concat_input,
            outdir
        )

        final_concat_primers_input = CONCAT_PRIMERS.out.concat_primers_out
        .groupTuple(by: [0, 1])


        FINAL_CONCAT_PRIMERS(
            final_concat_primers_input,
            outdir
        )

        // Join concatenated primers to the fastp-cleaned paired reads files and run cutadapt on them
        cutadapt_input = FINAL_CONCAT_PRIMERS.out.final_concat_primers_out
                        .join(fastp_cleaned_fastq, by: [0, 1])

        CUTADAPT(
            cutadapt_input,
            outdir
        )

    emit:
        final_concat_primers_out = FINAL_CONCAT_PRIMERS.out.final_concat_primers_out
        cutadapt_out = CUTADAPT.out.cutadapt_out
    
}