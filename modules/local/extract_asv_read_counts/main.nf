
process EXTRACT_ASV_READ_COUNTS {
    tag "$meta.id"
    label 'very_light'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "https://depot.galaxyproject.org/singularity/mgnify-pipelines-toolkit:${params.mpt_version}":
        "biocontainers/mgnify-pipelines-toolkit:${params.mpt_version}" }"
    
    input:
    tuple val(meta), path(silva_asvs), path(pr2_asvs)

    output:
    tuple val(meta), path("*asv_read_counts.tsv"), emit: asvs_left
    path "versions.yml"                          , emit: versions

    script:
    """
    cat ${silva_asvs} ${pr2_asvs} | sort | uniq > ${meta.id}_${meta.var_region}_asv_read_counts.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mgnify-pipelines-toolkit: ${params.mpt_version}
    END_VERSIONS
    """
}
