
process GENERAL_PRIMER_FLAG {
    // Check for the presence of primers in general
    tag "$meta.id"
    label 'very_light'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "https://depot.galaxyproject.org/singularity/mi-pimento:1.0.0--pyhdfd78af_0":
        "biocontainers/mi-pimento:1.0.0--pyhdfd78af_0"}"
    input:
    tuple val(meta), path(reads_merged)

    output:
    tuple val(meta), path("*general_primer_out.txt"), emit: general_primer_out
    path "versions.yml"                             , emit: versions

    script:
    """
    pimento are_there_primers -i ${reads_merged} -o ${meta.id}_${meta.var_region}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mi-pimento: \$( pimento --version | cut -d" " -f3 )
    END_VERSIONS
    """
}
