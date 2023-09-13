
include { CMSEARCH } from '../modules/cmsearch.nf'
include { CMSEARCH_DEOVERLAP } from '../modules/cmsearch_deoverlap.nf'
include { EASEL } from '../modules/easel.nf'
include { EXTRACT_COORDS } from '../modules/extract_coords.nf'

workflow CMSEARCH_SUBWF {

    // Subworkflow that runs cmsearch to find matches in RFAM in fasta files
    
    take:
        fasta
        outdir

    main:
        CMSEARCH(
            fasta, 
            outdir
        )

        CMSEARCH_DEOVERLAP(
            CMSEARCH.out.cmsearch_out,
            outdir
        )

        easel_input = fasta
                      .join(CMSEARCH_DEOVERLAP.out.cmsearch_deoverlap_out, by: [0, 1])
        
        EASEL(
            easel_input,
            outdir
        )

        EXTRACT_COORDS(
            EASEL.out.easel_coords,
            EASEL.out.matched_seqs_with_coords,
            outdir
        )

    emit:
        cmsearch_out = CMSEARCH.out.cmsearch_out
        cmsearch_deoverlap_out = CMSEARCH_DEOVERLAP.out.cmsearch_deoverlap_out
        easel_out = EASEL.out.easel_coords
        ssu_fasta = EXTRACT_COORDS.out.ssu_fasta
        lsu_fasta = EXTRACT_COORDS.out.lsu_fasta
        concat_ssu_lsu_coords = EXTRACT_COORDS.out.concat_ssu_lsu_coords
    
}