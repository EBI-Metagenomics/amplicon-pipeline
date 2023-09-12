
include { PRIMER_VALIDATION_SEARCH } from '../modules/primer_validation_search.nf'
include { PRIMER_VALIDATION_DEOVERLAP } from '../modules/primer_validation_deoverlap.nf'
include { PRIMER_VALIDATION_CLASSIFY_VAR_REGIONS } from '../modules/primer_validation_classify_var_regions.nf'

workflow PRIMER_VALIDATION {
    
    take:
        primer_validation_input
        outdir

    main:
        PRIMER_VALIDATION_SEARCH(
            primer_validation_input,
            outdir
        )

        PRIMER_VALIDATION_DEOVERLAP(
            PRIMER_VALIDATION_SEARCH.out.cmsearch_out,
            outdir
        )

        primer_validation_classify_var_regions_input = PRIMER_VALIDATION_DEOVERLAP.out.cmsearch_deoverlap_out
                                                    .join(primer_validation_input, by: [0, 1])

        PRIMER_VALIDATION_CLASSIFY_VAR_REGIONS(
            primer_validation_classify_var_regions_input,
            outdir
        )

    emit:
        primer_validation_out = PRIMER_VALIDATION_CLASSIFY_VAR_REGIONS.out.primer_validation_out
    
}