
include { FORMAT_BEDFILE } from '../../modules/local/format_bedfile/main.nf'
include { BEDTOOLS_MASKFASTA } from '../../modules/nf-core/bedtools/maskfasta/main'                                                                                                                          


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

        BEDTOOLS_MASKFASTA(
            bedtools_input
        )

    emit:
        its_masked_out = BEDTOOLS_MASKFASTA.out.fasta
    
}