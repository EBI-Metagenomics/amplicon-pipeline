
include { CLASSIFY_VAR_REGIONS } from '../modules/classify_var_regions.nf'
include { EXTRACT_VAR_REGIONS } from '../modules/extract_var_regions.nf'

workflow AMP_REGION_INFERENCE {
    
    take:
        cmsearch_deoverlap_out
        merged_reads
        outdir

    main:

        CLASSIFY_VAR_REGIONS(
            cmsearch_deoverlap_out,
            outdir
        )

        extract_var_input = CLASSIFY_VAR_REGIONS.out.classify_var_regions
        .transpose()
        .combine(merged_reads, by: [0, 1])

        EXTRACT_VAR_REGIONS(
            extract_var_input,
            outdir
        )

    emit:
        concat_var_regions = CLASSIFY_VAR_REGIONS.out.concat_var_regions
        extracted_var_out = EXTRACT_VAR_REGIONS.out.extracted_var_out
        extracted_var_path = EXTRACT_VAR_REGIONS.out.extracted_var_path
    
}