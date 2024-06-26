
process GENERAL_PRIMER_FLAG {
    // Check for the presence of primers in general
    tag "$meta.id"
    label 'very_light'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "https://depot.galaxyproject.org/singularity/mgnify-pipelines-toolkit:${params.mpt_version}":
        "biocontainers/mgnify-pipelines-toolkit:${params.mpt_version}" }"

    input:
    tuple val(meta), path(reads_merged)

    output:
    tuple val(meta), path("*general_primer_out.txt"), emit: general_primer_out
    path "versions.yml"                             , emit: versions

    script:
    """
    are_there_primers -i $reads_merged -s ${meta.id}_${meta.var_region} -o ./

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mgnify-pipelines-toolkit: ${params.mpt_version}
    END_VERSIONS
    """
}
