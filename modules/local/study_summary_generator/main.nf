
process STUDY_SUMMARY_GENERATOR {
    label 'very_light'
    container "community.wave.seqera.io/library/mgnify-pipelines-toolkit_pip_pyfastx:b75e1f542397fdb7" // TODO: temporary wave container until bugfix

    input:
    path pipeline_outputs
    path successful_runs_file
    val prefix
    val non_insdc

    output:
    path "*.tsv"        , emit: tax_study_summary
    path "*.json"       , emit: markergene_study_summaries, optional:true
    path "versions.yml" , emit: versions

    script:
    def insdc_input = non_insdc ? "--non_insdc" : ""

    """
    # TODO: fix incoming for this in the toolkit, can remove it soon
    python /opt/conda/lib/python3.11/site-packages/mgnify_pipelines_toolkit/analysis/shared/study_summary_generator.py summarise -a ${pipeline_outputs} -r ${successful_runs_file} -p ${prefix} ${insdc_input}

    markergene_study_summary.py -i ${pipeline_outputs} -r ${successful_runs_file} -p ${prefix}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mgnify-pipelines-toolkit: ${params.mpt_version}
    END_VERSIONS
    """

}
