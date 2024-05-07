
process GENERAL_PRIMER_FLAG {
    // Check for the presence of primers in general
    tag "$meta.id"
    label 'very_light'
    conda "${projectDir}/conf/environment.yml"

    input:
    tuple val(meta), path(reads_merged)

    output:
    tuple val(meta), path("*general_primer_out.txt"), emit: general_primer_out

    """
    are_there_primers -i $reads_merged -s ${meta.id}_${meta.var_region} -o ./
    """

}
