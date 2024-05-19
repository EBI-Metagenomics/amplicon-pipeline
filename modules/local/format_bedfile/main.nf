
process FORMAT_BEDFILE {
    tag "$meta.id"
    label 'very_light'

    // TODO: this one should include a container, there are several versions of awk
    // and not all of them have the same flags, this one should be fine but it's
    // a good practice
    input:
    tuple val(meta), path(concat_ssu_lsu_coords)

    output:
    tuple val(meta), path("*maskfile.bed"), emit: format_bedfile_out

    script:
    """
    awk '\$2 > \$3 { var = \$3; \$3 = \$2; \$2 = var } 1 {print \$4,\$2,\$3}' OFS='\t' $concat_ssu_lsu_coords > ${meta.id}_maskfile.bed
    """

}
