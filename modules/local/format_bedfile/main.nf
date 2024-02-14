
process FORMAT_BEDFILE {
    tag "$meta.id"
    label 'light'

    input:
    tuple val(meta), path(concat_ssu_lsu_coords)

    output:
    tuple val(meta), path("*maskfile.bed"), emit: format_bedfile_out

    """
    awk '\$2 > \$3 { var = \$3; \$3 = \$2; \$2 = var } 1 {print \$4,\$2,\$3}' OFS='\t' $concat_ssu_lsu_coords > ${meta.id}_maskfile.bed
    """

}
