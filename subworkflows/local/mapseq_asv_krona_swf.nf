
include { MAPSEQ                } from '../../modules/ebi-metagenomics/mapseq/main'
include { MAPSEQ2ASVTABLE       } from '../../modules/local/mapseq2asvtable/main.nf'
include { MAKE_ASV_COUNT_TABLES } from '../../modules/local/make_asv_count_tables/main.nf'
include { KRONA_KTIMPORTTEXT    } from '../../modules/ebi-metagenomics/krona/ktimporttext/main'
include { EXTRACT_ASVS_LEFT     } from '../../modules/local/extract_asvs_left/main.nf'

workflow MAPSEQ_ASV_KRONA {
    
    take:
        dada2_output
        concat_var_regions
        extracted_var_path
        krona_tuple

    main:

        ch_versions = channel.empty()

        mapseq_input = dada2_output
            .map { meta, _maps, asv_seqs, _filt_reads ->
                [ meta, asv_seqs ]
            }
            .combine (
                krona_tuple.map { fasta, tax, otu, _mscluster, _label -> 
                    [fasta, tax, otu] 
                }
            )
            .multiMap{ meta, input, reference ->
                input: [meta, input]
                reference: [meta, reference]
            }
        MAPSEQ(
            mapseq_input.input,
            mapseq_input.reference,
        )
        ch_versions = ch_versions.mix(MAPSEQ.out.versions.first())

        mapseq2asvtable_input = MAPSEQ.out.mseq
            .combine (
                krona_tuple.map { _fasta, _tax, _otu, _mscluster, label -> 
                    label 
                }
            )
            .multiMap{ meta, mseq, ref_label ->
                mseq: [meta, mseq]
                ref_label: [meta, ref_label]
            }
        MAPSEQ2ASVTABLE(
            mapseq2asvtable_input.mseq,
            mapseq2asvtable_input.ref_label
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
                      .map { meta, maps, _asv_seqs, filt_reads -> 
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
                               .map { meta, _var_region, var_regions_size, maps, asvtaxtable, filt_reads, __, concat_str, concat_vars ->
                                    [ meta + ['var_regions_size':var_regions_size], concat_str, maps, asvtaxtable, filt_reads, concat_vars ]
                               }
        // Add in the concatenated var region channel to the rest of the input
        final_asv_count_table_input = split_input
                                      .mix(multi_region_concats)
                                      .map { meta, var_region, maps, asvtaxtable, filt_reads, extracted_var ->
                                        [ meta + ['var_region': var_region], maps, asvtaxtable, filt_reads, extracted_var ]
                                      }
        
        make_table_input = final_asv_count_table_input
            .combine (
                krona_tuple.map { _fasta, _tax, _otu, _mscluster, label -> 
                    label 
                }
            )
            .multiMap{ meta, input, ref_label ->
                input: [meta, input]
                ref_label: [meta, ref_label]
            }
        MAKE_ASV_COUNT_TABLES(
            make_table_input.input,
            make_table_input.ref_label
        )
        ch_versions = ch_versions.mix(MAKE_ASV_COUNT_TABLES.out.versions.first())

        KRONA_KTIMPORTTEXT(
            MAKE_ASV_COUNT_TABLES.out.asv_krona_counts,
        )
        ch_versions = ch_versions.mix(KRONA_KTIMPORTTEXT.out.versions.first())

    emit:
        asv_krona_counts = MAKE_ASV_COUNT_TABLES.out.asv_krona_counts
        asv_read_counts = MAKE_ASV_COUNT_TABLES.out.asv_read_counts
        asvtaxtable = MAPSEQ2ASVTABLE.out.asvtaxtable
        krona_out = KRONA_KTIMPORTTEXT.out.html
        versions = ch_versions
}
