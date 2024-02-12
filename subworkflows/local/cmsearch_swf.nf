
include { RRNA_EXTRACTION } from '../../subworkflows/ebi-metagenomics/rrna_extraction/main'
include { EXTRACT_COORDS } from '../../modules/local/extract_coords/main.nf'

workflow CMSEARCH_SUBWF {

    // Subworkflow that runs cmsearch to find matches in RFAM in fasta files
    
    take:
        reads_fasta
        rfam
        claninfo
    main:

        RRNA_EXTRACTION(
            reads_fasta,
            rfam,
            claninfo
        )

        EXTRACT_COORDS(
            RRNA_EXTRACTION.out.easel_sfetch,
            RRNA_EXTRACTION.out.matched_seqs_with_coords,
        )

    emit:
        cmsearch_deoverlap_out = RRNA_EXTRACTION.out.cmsearch_deoverlap
        easel_out = RRNA_EXTRACTION.out.easel_sfetch
        ssu_fasta = EXTRACT_COORDS.out.ssu_fasta
        lsu_fasta = EXTRACT_COORDS.out.lsu_fasta
        concat_ssu_lsu_coords = EXTRACT_COORDS.out.concat_ssu_lsu_coords
    
}