
include { CONCAT_PRIMERS } from '../../modules/local/concat_primers.nf'
include { FINAL_CONCAT_PRIMERS } from '../../modules/local/final_concat_primers.nf'
include { REV_COMP_SE_PRIMERS } from '../../modules/local/rev_comp_se_primers.nf'
include { CUTADAPT } from '../../modules/local/cutadapt.nf'

workflow CONCAT_PRIMER_CUTADAPT {
    
    take:
        concat_input
        reads
    main:
        CONCAT_PRIMERS(
            concat_input
        )

        final_concat_primers_input = CONCAT_PRIMERS.out.concat_primers_out
        .groupTuple(by: 0)

        FINAL_CONCAT_PRIMERS(
            final_concat_primers_input
        )

        REV_COMP_SE_PRIMERS(
            FINAL_CONCAT_PRIMERS.out.final_concat_primers_out
        )

        // Join concatenated primers to the fastp-cleaned paired reads files and run cutadapt on them
        cutadapt_input = REV_COMP_SE_PRIMERS.out.rev_comp_se_primers_out
                        .join(reads, by: [0])

        CUTADAPT(
            cutadapt_input
        )

    emit:
        final_concat_primers_out = FINAL_CONCAT_PRIMERS.out.final_concat_primers_out
        cutadapt_out = CUTADAPT.out.cutadapt_out
    
}