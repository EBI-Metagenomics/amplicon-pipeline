
process STUDY_SUMMARY_GENERATOR {
    tag "$meta.id"
    label 'very_light'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "https://depot.galaxyproject.org/singularity/mgnify-pipelines-toolkit:${params.mpt_version}":
        "biocontainers/mgnify-pipelines-toolkit:${params.mpt_version}" }"

    input:
    path pipeline_outputs
    path successful_runs_file
    val prefix
    val non_insdc

    output:
    path "*.tsv"        , emit: study_summaries
    path "versions.yml" , emit: versions

    script:
    """
    study_summary_generator summarise -a ${pipeline_outputs} -r ${successful_runs_file} -p ${prefix} ${non_insdc}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mgnify-pipelines-toolkit: ${params.mpt_version}
    END_VERSIONS
    """

}
