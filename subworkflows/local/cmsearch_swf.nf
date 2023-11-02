
include { INFERNAL_CMSEARCH } from '../../modules/ebi-metagenomics/infernal/cmsearch/main.nf'
include { CMSEARCHTBLOUTDEOVERLAP } from '../../modules/ebi-metagenomics/cmsearchtbloutdeoverlap/main.nf'
include { EASEL } from '../../modules/local/easel/main.nf'
include { EXTRACT_COORDS } from '../../modules/local/extract_coords.nf'

workflow CMSEARCH_SUBWF {

    // Subworkflow that runs cmsearch to find matches in RFAM in fasta files
    
    take:
        reads_fasta
    main:

        INFERNAL_CMSEARCH(
            reads_fasta,
            file(params.rfam)
        )

        CMSEARCHTBLOUTDEOVERLAP(
            INFERNAL_CMSEARCH.out.cmsearch_tbl,
            file(params.rfam_clan)
        )
        
        ch_easel_input = reads_fasta
                         .join(CMSEARCHTBLOUTDEOVERLAP.out.cmsearch_tblout_deoverlapped)

        ch_easel_input.view()

    //     EASEL(
    //         ch_easel_input,
    //     )

    //     EXTRACT_COORDS(
    //         EASEL.out.easel_coords,
    //         EASEL.out.matched_seqs_with_coords,
    //     )

    // emit:
    //     cmsearch_out = INFERNAL_CMSEARCH.out.cmsearch_tbl
    //     cmsearch_deoverlap_out = CMSEARCHTBLOUTDEOVERLAP.out.cmsearch_tblout_deoverlapped
    //     easel_out = EASEL.out.easel_coords
    //     ssu_fasta = EXTRACT_COORDS.out.ssu_fasta
    //     lsu_fasta = EXTRACT_COORDS.out.lsu_fasta
    //     concat_ssu_lsu_coords = EXTRACT_COORDS.out.concat_ssu_lsu_coords
    
}