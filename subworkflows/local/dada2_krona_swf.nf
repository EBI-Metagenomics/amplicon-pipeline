
include { REMOVE_AMBIGUOUS_READS } from '../../modules/local/remove_ambiguous_reads/main.nf'
include { DADA2 } from '../../modules/local/dada2/main.nf'
include { MAKE_ASV_COUNT_TABLES } from '../../modules/local/make_asv_count_tables/main.nf'
include { KRONA } from '../../modules/local/krona.nf'

workflow DADA2_KRONA {
    
    take:
        dada2_input
        concat_var_regions
        extracted_var_path
        fastp_cleaned_fastq
        dada2_db
        krona_tuple

    main:

        REMOVE_AMBIGUOUS_READS(
            dada2_input
        )

        DADA2(
            REMOVE_AMBIGUOUS_READS.out.noambig_out,
            dada2_db,
            krona_tuple[4] // db_label
        )

        split_input = DADA2.out.dada2_out
                      .map { meta, maps, taxa -> 
                        [ meta.subMap('id', 'single_end', 'var_regions_size'), meta['var_region'], maps, taxa ]
                       }
                      .transpose(by: 1)
                      .join(extracted_var_path, by: [0, 1])
        
        multi_region_concats = split_input
                                .map { meta, var_region, maps, taxa, extracted_var ->
                                    [ meta.subMap('id', 'single_end'), var_region, meta['var_regions_size'], maps, taxa, extracted_var ]
                                }
                               .join(concat_var_regions, by: 0)
                               .map {meta, var_region, var_regions_size, maps, taxa, _, concat_str, concat_vars ->
                                    [ meta + ['var_regions_size':var_regions_size], concat_str, maps, taxa, concat_vars ]
                               }


        final_asv_count_table_input = split_input
                                      .mix(multi_region_concats)
                                      .map { meta, var_region, maps, taxa, extracted_var ->
                                        [ meta.subMap('id', 'single_end'), var_region, meta['var_regions_size'], maps, taxa, extracted_var ]
                                    }
                                      .combine(fastp_cleaned_fastq, by: 0)
                                      .map { meta, var_region, var_regions_size, maps, taxa, extracted_var, reads ->
                                        [ meta + ['var_region': var_region, 'var_regions_size': var_regions_size], maps, taxa, extracted_var, reads ]
                                      }    
        MAKE_ASV_COUNT_TABLES(
            final_asv_count_table_input
        )

        KRONA(
            MAKE_ASV_COUNT_TABLES.out.asv_count_tables_out,
            krona_tuple
        )

    emit:
        asv_count_tables_out = MAKE_ASV_COUNT_TABLES.out.asv_count_tables_out
        krona_out = KRONA.out.krona_out
    
}