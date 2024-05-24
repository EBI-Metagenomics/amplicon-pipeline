
include { FORMAT_BEDFILE     } from '../../modules/local/format_bedfile/main.nf'
include { BEDTOOLS_MASKFASTA } from '../../modules/nf-core/bedtools/maskfasta/main'


workflow MASK_FASTA_SWF {
    
    take:
        reads_fasta
        concat_ssu_lsu_coords
    main:

        ch_versions = Channel.empty()

        FORMAT_BEDFILE(
            concat_ssu_lsu_coords,
        )
        ch_versions = ch_versions.mix(FORMAT_BEDFILE.out.versions.first())

        bedtools_input = reads_fasta
                         .join(FORMAT_BEDFILE.out.format_bedfile_out)

        BEDTOOLS_MASKFASTA(
            bedtools_input
        )
        ch_versions = ch_versions.mix(BEDTOOLS_MASKFASTA.out.versions.first())

    emit:
        masked_out = BEDTOOLS_MASKFASTA.out.fasta
        versions = ch_versions
}