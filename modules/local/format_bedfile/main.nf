
process FORMAT_BEDFILE {
    tag "$meta.id"
    label 'very_light'
    container 'docker://quay.io/biocontainers/mgnify-pipelines-toolkit:0.1.2--pyhdfd78af_0'

    input:
    tuple val(meta), path(concat_ssu_lsu_coords)

    output:
    tuple val(meta), path("*maskfile.bed"), emit: format_bedfile_out

    script:
    """
    awk '\$2 > \$3 { var = \$3; \$3 = \$2; \$2 = var } 1 {print \$4,\$2,\$3}' OFS='\t' $concat_ssu_lsu_coords > ${meta.id}_maskfile.bed
    """

}
