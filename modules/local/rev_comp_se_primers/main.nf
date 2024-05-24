
process REV_COMP_SE_PRIMERS {
    tag "$meta.id"
    label 'very_light'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "https://depot.galaxyproject.org/singularity/mgnify-pipelines-toolkit:${params.mpt_version}":
        "biocontainers/mgnify-pipelines-toolkit:${params.mpt_version}" }"

    input:
    tuple val(meta), path(final_concat_primers)

    output:
    tuple val(meta), path("*rev_comp_se_primers.fasta"), optional: true, emit: rev_comp_se_primers_out
    path "versions.yml"                                , emit: versions

    script:
    """
    if [[ -s $final_concat_primers && ${meta.single_end} = true ]]; then
        rev_comp_se_primers -i $final_concat_primers -s ${meta.id} -o ./
    elif [[ -s $final_concat_primers && ${meta.single_end} = false ]]; then
        cat $final_concat_primers > ./${meta.id}_rev_comp_se_primers.fasta
    else
        touch ./${meta.id}_rev_comp_se_primers.fasta
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mgnify-pipelines-toolkit: ${params.mpt_version}
    END_VERSIONS
    """

}
