
include { PRIMER_VALIDATION_SEARCH } from '../../modules/local/primer_validation_search.nf'
include { PRIMER_VALIDATION_DEOVERLAP } from '../../modules/local/primer_validation_deoverlap.nf'
include { PRIMER_VALIDATION_CLASSIFY_VAR_REGIONS } from '../../modules/local/primer_validation_classify_var_regions/main.nf'

workflow PRIMER_VALIDATION {
    
    take:
        primer_validation_input
    main:
        PRIMER_VALIDATION_SEARCH(
            primer_validation_input,
            file(params.rfam)
        )

        PRIMER_VALIDATION_DEOVERLAP(
            PRIMER_VALIDATION_SEARCH.out.cmsearch_out,
            file(params.rfam_clan)
        )

        primer_validation_classify_var_regions_input = PRIMER_VALIDATION_DEOVERLAP.out.cmsearch_deoverlap_out
                                                       .join(primer_validation_input, by: 0)

        primer_validation_classify_var_regions_input.view()
        
        PRIMER_VALIDATION_CLASSIFY_VAR_REGIONS(
            primer_validation_classify_var_regions_input
        )

    emit:
        primer_validation_out = PRIMER_VALIDATION_CLASSIFY_VAR_REGIONS.out.primer_validation_out
    
}