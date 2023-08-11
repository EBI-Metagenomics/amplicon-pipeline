
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
                      .join(CMSEARCH_DEOVERLAP.out.cmsearch_deoverlap_out)
        
        EASEL(
            easel_input,
            outdir
        )

        EXTRACT_COORDS(
            EASEL.out.easel_coords,
            outdir
        )

    emit:
        cmsearch_out = CMSEARCH.out.cmsearch_out
        cmsearch_deoverlap_out = CMSEARCH_DEOVERLAP.out.cmsearch_deoverlap_out
        easel_out = EASEL.out.easel_coords
        ssu_fasta = EXTRACT_COORDS.out.ssu_fasta
        lsu_fasta = EXTRACT_COORDS.out.lsu_fasta
    
}