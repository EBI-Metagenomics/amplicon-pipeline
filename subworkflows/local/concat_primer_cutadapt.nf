
include { CONCAT_PRIMERS } from '../../modules/local/concat_primers/main.nf'
include { FINAL_CONCAT_PRIMERS } from '../../modules/local/final_concat_primers/main.nf'
include { REV_COMP_SE_PRIMERS } from '../../modules/local/rev_comp_se_primers/main.nf'
include { SPLIT_PRIMERS_BY_STRAND } from '../../modules/local/split_primers_by_strand.nf'
include { CUTADAPT } from '../../modules/ebi-metagenomics/cutadapt/main.nf'

workflow CONCAT_PRIMER_CUTADAPT {
    
    take:
        concat_input
        reads
    main:

        CONCAT_PRIMERS(
            concat_input
        )

        final_concat_primers_input = CONCAT_PRIMERS.out.concat_primers_out
                                    .map{ tuple(["id":it[0].id, "single_end":it[0].single_end], it[0].var_region, it[1]) }
                                    .groupTuple(by: 0)
                                    .map { tuple(it[0] + ["var_region":it[1]], it[2]) }

        FINAL_CONCAT_PRIMERS(
            final_concat_primers_input
        )

        REV_COMP_SE_PRIMERS(
            FINAL_CONCAT_PRIMERS.out.final_concat_primers_out
        )

        SPLIT_PRIMERS_BY_STRAND(
            REV_COMP_SE_PRIMERS.out.rev_comp_se_primers_out
        )

        // Join concatenated primers to the fastp-cleaned paired reads files and run cutadapt on them
        cutadapt_input = SPLIT_PRIMERS_BY_STRAND.out.stranded_primer_out
                        .map{ tuple(["id":it[0].id, "single_end":it[0].single_end], it[0].var_region, [it[1], it[2]]) }
                        .join(reads, by: [0])
                        .map{ tuple(it[0], it[3], it[2]) }

        CUTADAPT(
            cutadapt_input
        )

    emit:
        final_concat_primers_out = FINAL_CONCAT_PRIMERS.out.final_concat_primers_out
        cutadapt_out = CUTADAPT.out.reads
    
}