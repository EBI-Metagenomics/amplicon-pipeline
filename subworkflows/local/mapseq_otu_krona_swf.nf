
include { MAPSEQ } from '../../modules/local/mapseq.nf'
include { MAPSEQ2BIOM } from '../../modules/local/mapseq2biom.nf'
include { KRONA } from '../../modules/local/krona.nf'

workflow MAPSEQ_OTU_KRONA {
    
    take:
        subunit_fasta
        db_tuple
    main:

        MAPSEQ(
            subunit_fasta,
            db_tuple
        )

        MAPSEQ2BIOM(
            MAPSEQ.out.mapseq_out,
            db_tuple
        )

        KRONA(
            MAPSEQ2BIOM.out.mapseq2biom_krona_out,
            db_tuple
        )

    emit:
        krona_input = MAPSEQ2BIOM.out.mapseq2biom_krona_out
        krona_out = KRONA.out.krona_out
    
}