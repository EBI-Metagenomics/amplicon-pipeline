
process PERMUTE_PRIMERS {
    tag "$meta.id"
    label 'very_light'
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     "https://depot.galaxyproject.org/singularity/mgnify-pipelines-toolkit:${params.mpt_version}":
    //     "biocontainers/mgnify-pipelines-toolkit:${params.mpt_version}" }"
    container "oras://community.wave.seqera.io/library/pip_mgnify-pipelines-toolkit:64432eed2c673687"

    input:
    tuple val(meta), path(concat_primers_fasta)

    output:
    tuple val(meta), path("*permuted_primers.fasta"), emit: permuted_primers
    path "versions.yml"                             , emit: versions

    script:
    """
    permute_primers -i $concat_primers_fasta -p ${meta.id}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mgnify-pipelines-toolkit: ${params.mpt_version}
    END_VERSIONS
    """
}
