
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
                      .map{ tuple(["id":it[0].id, "single_end":it[0].single_end], it[0].var_region, it[1], it[2], it[3], it[4]) }
                      .transpose(by: 1)
                      .join(extracted_var_path, by: [0, 1])

        multi_region_concats = split_input
                               .join(concat_var_regions, by: 0)
                               .map( {tuple(it[0], "concat", it[2], it[3], it[4], it[5], it[8])} )

        final_asv_count_table_input = split_input
                                      .mix(multi_region_concats)
                                      .combine(fastp_cleaned_fastq, by: 0)
                                      .map { tuple(it[0] + ["var_region":it[1]], it[2], it[3], it[4], it[5], it[6], it[7]) }        
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