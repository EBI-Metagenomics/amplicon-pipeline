
process FORMAT_BEDFILE {

    label 'light'

    input:
    tuple val(project), val(sampleId), path(concat_ssu_lsu_coords)
    val outdir

    output:
    tuple val(project), val(sampleId), path("*maskfile.bed"), emit: format_bedfile_out

    """
    awk '\$2 > \$3 { var = \$3; \$3 = \$2; \$2 = var } 1 {print \$4,\$2,\$3}' OFS='\t' $concat_ssu_lsu_coords > ${sampleId}_maskfile.bed
    """

}
