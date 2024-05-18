
include { REMOVE_AMBIGUOUS_READS } from '../../modules/local/remove_ambiguous_reads/main.nf'
include { DADA2                  } from '../../modules/local/dada2/main.nf'

workflow DADA2_SWF {
    
    take:
        dada2_input

    main:

        REMOVE_AMBIGUOUS_READS(
            dada2_input
        )

        DADA2(
            REMOVE_AMBIGUOUS_READS.out.noambig_out
        )

    emit:
        dada2_out = DADA2.out.dada2_out
    
}