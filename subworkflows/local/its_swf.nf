
include { FORMAT_BEDFILE } from '../../modules/local/format_bedfile.nf'
include { BEDTOOLS } from '../../modules/local/bedtools.nf'

workflow ITS_SWF {
    
    take:
        reads_fasta
        concat_ssu_lsu_coords
    main:

        FORMAT_BEDFILE(
            concat_ssu_lsu_coords,
        )

        bedtools_input = reads_fasta
                         .join(FORMAT_BEDFILE.out.format_bedfile_out)

        BEDTOOLS(
            bedtools_input,
        )

    emit:
        its_masked_out = BEDTOOLS.out.its_masked_out
    
}