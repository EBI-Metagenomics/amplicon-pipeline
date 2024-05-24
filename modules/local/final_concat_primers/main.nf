
process FINAL_CONCAT_PRIMERS {
    tag "$meta.id"
    label 'very_light'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "https://depot.galaxyproject.org/singularity/mgnify-pipelines-toolkit:${params.mpt_version}":
        "biocontainers/mgnify-pipelines-toolkit:${params.mpt_version}" }"
    
    input:
    tuple val(meta), path(cat_primers)

    output:
    tuple val(meta), path("*final_concat_primers.fasta"), optional: true, emit: final_concat_primers_out
    path "versions.yml"                                 , emit: versions

    script:
    """
    cat *concat_primers.fasta > temp_concat_primers.fasta
    if [[ -s temp_concat_primers.fasta ]]; then
        awk '!a[\$0]++' temp_concat_primers.fasta > ${meta.id}_final_concat_primers.fasta
    else
       touch ${meta.id}_final_concat_primers.fasta
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mgnify-pipelines-toolkit: ${params.mpt_version}
    END_VERSIONS
    """

}
