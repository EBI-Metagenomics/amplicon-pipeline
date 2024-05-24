
process CONCAT_PRIMERS {
    // Concatenate result of both the std primer identification and auto primer identification
    // This is important in case for a run the fwd strand has a known std primer but the reverse
    // strand has an auto primer
    tag "$meta.id"
    label 'extremely_light'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "https://depot.galaxyproject.org/singularity/mgnify-pipelines-toolkit:${params.mpt_version}":
        "biocontainers/mgnify-pipelines-toolkit:${params.mpt_version}" }"

    input:
    tuple val(meta), path(std_primers), path(auto_primers)

    output:
    tuple val(meta), path("*concat_primers.fasta") , optional: true, emit: concat_primers_out
    path "versions.yml"                            , emit: versions

    script:
    """
    if [[ -s ${auto_primers} ]]
    then
        cat ${std_primers} > ${meta.id}_${meta.var_region}_concat_primers.fasta
        echo -en "\\n" >> ${meta.id}_${meta.var_region}_concat_primers.fasta
        cat ${auto_primers} >> ${meta.id}_${meta.var_region}_concat_primers.fasta
    else
        cat ${std_primers} > ${meta.id}_${meta.var_region}_concat_primers.fasta
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mgnify-pipelines-toolkit: ${params.mpt_version}
    END_VERSIONS
    """

}
