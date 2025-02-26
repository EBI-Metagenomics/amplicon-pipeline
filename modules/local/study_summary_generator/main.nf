
process STUDY_SUMMARY_GENERATOR {
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
    def insdc_input = non_insdc ? "--non_insdc" : ""

    """
    # TODO: fix incoming for this in the toolkit, can remove it soon
    python /usr/local/lib/python3.11/site-packages/mgnify_pipelines_toolkit/analysis/shared/study_summary_generator.py summarise -a ${pipeline_outputs} -r ${successful_runs_file} -p ${prefix} ${insdc_input}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mgnify-pipelines-toolkit: ${params.mpt_version}
    END_VERSIONS
    """

}
