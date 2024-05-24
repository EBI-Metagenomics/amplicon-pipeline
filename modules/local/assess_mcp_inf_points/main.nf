
process ASSESS_MCP_INF_POINTS {
    // Select inflection points most likely to be primer cutoff points
    tag "$meta.id"
    label 'very_light'
    // TODO: uncomment container when you release fix to mpt
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     "https://depot.galaxyproject.org/singularity/mgnify-pipelines-toolkit:${params.mpt_version}":
    //     "biocontainers/mgnify-pipelines-toolkit:${params.mpt_version}" }"

    input:
    tuple val(meta), path(inf_points_out), path(reads_merged)

    output:
    tuple val(meta), path("*auto_primers.fasta")    , emit: auto_primer_out
    path "versions.yml"                             , emit: versions

    script:
    """
    if [[ -s ./$inf_points_out ]]; then
        assess_inflection_point_mcp -i $reads_merged -p $inf_points_out -s ${meta.id}_${meta.var_region} -o ./
    else
        touch ${meta.id}_${meta.var_region}_auto_primers.fasta
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mgnify-pipelines-toolkit: ${params.mpt_version}
    END_VERSIONS
    """

}
