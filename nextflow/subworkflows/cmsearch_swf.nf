
include { cmsearch } from '../modules/cmsearch.nf'

workflow CMSEARCH {

    take:
        fasta
        outdir

    main:
        cmsearch(fasta, outdir)
        // cmsearch_deoverlap
    
    emit:
        cmsearch_out = cmsearch.out.cmsearch_out
    
}