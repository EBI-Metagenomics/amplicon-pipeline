
include { FORMAT_BEDFILE } from '../modules/format_bedfile.nf'
include { BEDTOOLS } from '../modules/bedtools.nf'

workflow ITS_SWF {
    
    take:
        fasta
        concat_ssu_lsu_coords
        outdir

    main:

        FORMAT_BEDFILE(
            concat_ssu_lsu_coords,
            outdir
        )

        bedtools_input = fasta
                         .join(FORMAT_BEDFILE.out.format_bedfile_out, by: [0, 1])


        BEDTOOLS(
            bedtools_input,
            outdir
        )

    emit:
        its_masked_out = BEDTOOLS.out.its_masked_out
    
}