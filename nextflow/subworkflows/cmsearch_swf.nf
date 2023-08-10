
include { cmsearch } from '../modules/cmsearch.nf'
include { cmsearch_deoverlap } from '../modules/cmsearch_deoverlap.nf'
include { easel } from '../modules/easel.nf'

workflow CMSEARCH {

    // Subworkflow that runs cmsearch to find matches in RFAM in fasta files
    
    take:
        fasta
        outdir

    main:
        cmsearch(fasta, outdir)
        cmsearch_deoverlap(cmsearch.out.cmsearch_out, outdir)

        easel_input = fasta.join(cmsearch_deoverlap.out.cmsearch_deoverlap_out)
        easel(easel_input, outdir)
    
    emit:
        cmsearch_out = cmsearch.out.cmsearch_out
        cmsearch_deoverlap_out = cmsearch_deoverlap.out.cmsearch_deoverlap_out
        easel_out = easel.out.easel_coords
    
}