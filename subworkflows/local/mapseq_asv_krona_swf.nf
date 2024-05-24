
include { MAPSEQ                } from '../../modules/ebi-metagenomics/mapseq/main'
include { MAPSEQ2ASVTABLE       } from '../../modules/local/mapseq2asvtable/main.nf'
include { MAKE_ASV_COUNT_TABLES } from '../../modules/local/make_asv_count_tables/main.nf'
include { KRONA_KTIMPORTTEXT    } from '../../modules/ebi-metagenomics/krona/ktimporttext/main'


workflow MAPSEQ_ASV_KRONA {
    
    take:
        dada2_output
        concat_var_regions
        extracted_var_path
        krona_tuple

    main:

        ch_versions = Channel.empty()

        mapseq_input = dada2_output
                       .map { meta, maps, asv_seqs, filt_reads ->
                            [ meta, asv_seqs ]
                        }
        MAPSEQ(
            mapseq_input,
            tuple(krona_tuple[0], krona_tuple[1], krona_tuple[3])
        )
        ch_versions = ch_versions.mix(MAPSEQ.out.versions.first())

        MAPSEQ2ASVTABLE(
            MAPSEQ.out.mseq,
            krona_tuple[4] // db_label
        )
        ch_versions = ch_versions.mix(MAPSEQ2ASVTABLE.out.versions.first())

        // Transpose by var region in case any samples have more than one
        split_mapseq2asvtable = MAPSEQ2ASVTABLE.out.asvtaxtable
                                .map { meta, asvtaxtable -> 
                                     [ meta.subMap('id', 'single_end', 'var_regions_size'), meta['var_region'], asvtaxtable ]       
                                    }
                                .transpose(by: 1)

        // Transpose by var region in case any samples have more than one. Also reorder the inputs slightly
        split_input = dada2_output
                      .map { meta, maps, asv_seqs, filt_reads -> 
                        [ meta.subMap('id', 'single_end', 'var_regions_size'), meta['var_region'], maps, filt_reads ]
                       }
                      .transpose(by: 1)
                      .join(extracted_var_path, by: [0, 1])
                      .join(split_mapseq2asvtable, by: [0, 1])
                      .map {
                        meta, var_region, maps, filt_reads, extracted_var, asvtaxtable ->
                        [ meta, var_region, maps, asvtaxtable, filt_reads, extracted_var ] 
                      }

        // Make a channel containing the concatenated var region for any sample that has more than one var region
        multi_region_concats = split_input
                                .map { meta, var_region, maps, asvtaxtable, filt_reads, extracted_var ->
                                    [ meta.subMap('id', 'single_end'), var_region, meta['var_regions_size'], maps, asvtaxtable, filt_reads, extracted_var ]
                                }
                               .join(concat_var_regions, by: 0)
                               .map { meta, var_region, var_regions_size, maps, asvtaxtable, filt_reads, _, concat_str, concat_vars ->
                                    [ meta + ['var_regions_size':var_regions_size], concat_str, maps, asvtaxtable, filt_reads, concat_vars ]
                               }
        // Add in the concatenated var region channel to the rest of the input
        final_asv_count_table_input = split_input
                                      .mix(multi_region_concats)
                                      .map { meta, var_region, maps, asvtaxtable, filt_reads, extracted_var ->
                                        [ meta + ['var_region': var_region], maps, asvtaxtable, filt_reads, extracted_var ]
                                      }
        
        MAKE_ASV_COUNT_TABLES(
            final_asv_count_table_input
        )
        ch_versions = ch_versions.mix(MAKE_ASV_COUNT_TABLES.out.versions.first())

        KRONA_KTIMPORTTEXT(
            MAKE_ASV_COUNT_TABLES.out.asv_count_tables_out,
        )
        ch_versions = ch_versions.mix(KRONA_KTIMPORTTEXT.out.versions.first())

    emit:
        asv_count_tables_out = MAKE_ASV_COUNT_TABLES.out.asv_count_tables_out
        krona_out = KRONA_KTIMPORTTEXT.out.html
        versions = ch_versions
    
}