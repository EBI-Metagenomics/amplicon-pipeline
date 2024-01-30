
include { CLASSIFY_VAR_REGIONS } from '../../modules/local/classify_var_regions/main.nf'
include { SEQTK_SUBSEQ as EXTRACT_VAR_REGIONS } from '../../modules/nf-core/seqtk/subseq/main.nf'

workflow AMP_REGION_INFERENCE {
    
    take:
        cmsearch_deoverlap_out
        reads_merged
    main:

        CLASSIFY_VAR_REGIONS(
            cmsearch_deoverlap_out,
        )

        extract_var_input = CLASSIFY_VAR_REGIONS.out.classify_var_regions
                            .transpose()
                            .combine(reads_merged, by: 0)

        EXTRACT_VAR_REGIONS(
            extract_var_input,
        )

        final_var_out = EXTRACT_VAR_REGIONS.out.extracted_var_out
                    .map { tuple(it[0] + ["var_region":it[1]], it[2]) }
    emit:
        concat_var_regions = CLASSIFY_VAR_REGIONS.out.concat_var_regions
        extracted_var_out = final_var_out
        extracted_var_path = EXTRACT_VAR_REGIONS.out.extracted_var_path
    
}