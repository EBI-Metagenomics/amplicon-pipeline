
include { cmsearch } from '../modules/cmsearch.nf'
include { cmsearch_deoverlap } from '../modules/cmsearch_deoverlap.nf'
include { easel } from '../modules/easel.nf'
include { extract_coords } from '../modules/extract_coords.nf'

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
        extract_coords(easel.out.easel_coords, outdir)

    emit:
        cmsearch_out = cmsearch.out.cmsearch_out
        cmsearch_deoverlap_out = cmsearch_deoverlap.out.cmsearch_deoverlap_out
        easel_out = easel.out.easel_coords
        ssu_fasta = extract_coords.out.ssu_fasta
        lsu_fasta = extract_coords.out.lsu_fasta
    
}