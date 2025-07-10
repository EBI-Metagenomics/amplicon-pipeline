
include { REMOVE_AMBIGUOUS_READS } from '../../modules/local/remove_ambiguous_reads/main.nf'
include { DADA2                  } from '../../modules/local/dada2/main.nf'

workflow DADA2_SWF {
    
    take:
        dada2_input

    main:

        ch_versions = Channel.empty()

        REMOVE_AMBIGUOUS_READS(
            dada2_input
        )
        ch_versions = ch_versions.mix(REMOVE_AMBIGUOUS_READS.out.versions.first())

        DADA2(
            REMOVE_AMBIGUOUS_READS.out.noambig_out
        )
        ch_versions = ch_versions.mix(DADA2.out.versions.first())

    emit:
        dada2_out        = DADA2.out.dada2_out
        dada2_report     = DADA2.out.dada2_stats
        dada2_stats_fail = DADA2.out.dada2_stats_fail
        versions         = ch_versions
    
}