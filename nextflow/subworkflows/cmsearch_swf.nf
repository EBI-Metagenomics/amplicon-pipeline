
include { cmsearch } from '../modules/cmsearch.nf'

workflow CMSEARCH {

    // Subworkflow that runs cmsearch to find matches in RFAM in fasta files
    
    take:
        fasta
        outdir

    main:
        cmsearch(fasta, outdir)
        // cmsearch_deoverlap
    
    emit:
        cmsearch_out = cmsearch.out.cmsearch_out
    
}