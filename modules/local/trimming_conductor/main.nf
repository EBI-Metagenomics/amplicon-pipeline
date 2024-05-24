
process TRIMMING_CONDUCTOR {
    // Parse the outputs of the std and general primer searches and output
    // flags into "trimming_conductor_out.txt"
    tag "$meta.id"
    label 'very_light'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "https://depot.galaxyproject.org/singularity/mgnify-pipelines-toolkit:${params.mpt_version}":
        "biocontainers/mgnify-pipelines-toolkit:${params.mpt_version}" }"

    input:
    tuple val(meta), path(general_primer_flag), path(std_primer_flag)

    output:
    tuple val(meta), path("*trimming_conductor_out*"), emit: trimming_conductor_out
    path "versions.yml"                              , emit: versions

    script:
    """
    gen_fwd_flag="\$(sed '1q;d' ${general_primer_flag})"
    gen_rev_flag="\$(sed '2q;d' ${general_primer_flag})"

    std_fwd_primer="\$(grep 'F' ${std_primer_flag} || true)"
    std_rev_primer="\$(grep 'R' ${std_primer_flag} || true)"

    fwd_trim_flag='none'
    rev_trim_flag='none'

    if [[ ! -z \$std_fwd_primer ]]; then
        std_fwd_flag='1'
        fwd_trim_flag='std'
    else
        std_fwd_flag='0'
        if [[ \$gen_fwd_flag -eq '1' ]]; then
            fwd_trim_flag='auto'
        else
            fwd_trim_flag='none'
        fi
    fi

    if [[ ! -z \$std_rev_primer ]]; then
        std_rev_flag='1'
        rev_trim_flag='std'
    else
        std_rev_flag='0'
        if [[ \$gen_rev_flag -eq '1' ]]; then
            rev_trim_flag='auto'
        else
            rev_trim_flag='none'
        fi
    fi

    echo -en "\${fwd_trim_flag}\\n\${rev_trim_flag}" > ${meta.id}_trimming_conductor_out_${meta.var_region}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mgnify-pipelines-toolkit: ${params.mpt_version}
    END_VERSIONS
    """
}
