
include { CLASSIFY_VAR_REGIONS                } from '../../modules/local/classify_var_regions/main.nf'
include { SEQTK_SUBSEQ as EXTRACT_VAR_REGIONS } from '../../modules/nf-core/seqtk/subseq/main.nf'

workflow AMP_REGION_INFERENCE {
    
    take:
        cmsearch_deoverlap_out
        reads_merged
    main:

        ch_versions = channel.empty()

        CLASSIFY_VAR_REGIONS(
            cmsearch_deoverlap_out,
        )
        ch_versions = ch_versions.mix(CLASSIFY_VAR_REGIONS.out.versions.first())

        extract_var_input = CLASSIFY_VAR_REGIONS.out.classify_var_regions
                            .map { meta, var_regions -> 
                                [ meta, var_regions, var_regions.size() ]
                            }
                            .transpose()
                            .combine(reads_merged, by: 0)
                            .map { meta, var_regions, var_regions_size, reads ->
                                [ meta + ["var_regions_size":var_regions_size], var_regions, reads ]
                            }

        EXTRACT_VAR_REGIONS(
            extract_var_input,
        )
        ch_versions = ch_versions.mix(EXTRACT_VAR_REGIONS.out.versions.first())

        final_var_out = EXTRACT_VAR_REGIONS.out.extracted_var_out
                    .map { meta, var_region, extracted_reads ->
                        [ meta + ["var_region":var_region], extracted_reads]
                    }

    emit:
        concat_var_regions = CLASSIFY_VAR_REGIONS.out.concat_var_regions
        extracted_var_out = final_var_out
        extracted_var_path = EXTRACT_VAR_REGIONS.out.extracted_var_path
        versions = ch_versions
    
}
