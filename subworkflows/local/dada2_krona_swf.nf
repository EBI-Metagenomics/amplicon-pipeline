
include { REMOVE_AMBIGUOUS_READS } from '../modules/remove_ambiguous_reads.nf'
include { DADA2 } from '../modules/dada2.nf'
include { MAKE_ASV_COUNT_TABLES } from '../modules/make_asv_count_tables.nf'
include { KRONA } from '../modules/krona.nf'

workflow DADA2_KRONA {
    
    take:
        dada2_input
        concat_var_regions
        extracted_var_path
        fastp_cleaned_fastq
        silva_dada2_db
        ssu_mapseq_krona_tuple
        outdir

    main:
        REMOVE_AMBIGUOUS_READS(
            dada2_input,
            outdir
        )

        DADA2(
            REMOVE_AMBIGUOUS_READS.out.noambig_out,
            silva_dada2_db,
            outdir
        )

        split_input = DADA2.out.dada2_out
                      .transpose()
                      .join(extracted_var_path, by: [0, 1, 2])
                    

        multi_region_concats = split_input
                               .join(concat_var_regions, by: [0, 1])
                               .map( {tuple(it[0], it[1], "concat", it[3], it[4], it[5], it[6], it[7], it[10])} )
        
        final_asv_count_table_input = split_input
                                      .mix(multi_region_concats)
                                      .combine(fastp_cleaned_fastq, by: [0, 1])
                                    
        MAKE_ASV_COUNT_TABLES(
            final_asv_count_table_input,
            outdir
        )

        asv_krona_input = MAKE_ASV_COUNT_TABLES.out.asv_count_tables_out
                          .map( {it[0, 1, 3]} )
        KRONA(
            asv_krona_input,
            ssu_mapseq_krona_tuple,
            outdir
        )

    emit:
        asv_krona_input = asv_krona_input
        krona_out = KRONA.out.krona_out
    
}