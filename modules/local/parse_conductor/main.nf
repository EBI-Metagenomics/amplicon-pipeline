
process PARSE_CONDUCTOR {
    
    // Parse the trimming_conductor output and store 
    // the flags into environment variables 
    tag "$meta.id"
    label 'very_light'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "https://depot.galaxyproject.org/singularity/mgnify-pipelines-toolkit:${params.mpt_version}":
        "biocontainers/mgnify-pipelines-toolkit:${params.mpt_version}" }"

    input:
    tuple val(meta), path(trimming_conductor_out)

    output:
    tuple val(meta), env(fwd_flag), env(rev_flag), emit: conductor_out

    script:
    """
    fwd_flag=\$(sed '1q;d' "${trimming_conductor_out}")
    rev_flag=\$(sed '2q;d' "${trimming_conductor_out}")
    """
}
