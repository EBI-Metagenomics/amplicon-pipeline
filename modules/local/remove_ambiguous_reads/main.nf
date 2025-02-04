
process REMOVE_AMBIGUOUS_READS {
    // Run DADA2 pipeline including read-tracking
    tag "$meta.id"
    label 'process_medium'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "https://depot.galaxyproject.org/singularity/mgnify-pipelines-toolkit:${params.mpt_version}":
        "biocontainers/mgnify-pipelines-toolkit:${params.mpt_version}" }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*noambig*fastq.gz"), emit: noambig_out
    path "versions.yml"                       , emit: versions

    script:
    """
    if [[ ${meta.single_end} = true ]]; then
        remove_ambiguous_reads -f $reads -s ${meta.id}
    else 
        remove_ambiguous_reads -f ${reads[0]} -r ${reads[1]} -s ${meta.id}
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mgnify-pipelines-toolkit: ${params.mpt_version}
    END_VERSIONS
    """
}
