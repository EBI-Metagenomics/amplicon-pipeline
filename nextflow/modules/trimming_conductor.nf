
process trimming_conductor {
    // Parse the outputs of the std and general primer searches and output
    // flags into "trimming_conductor_out.txt"

    label 'light'
    publishDir "${outdir}/${project}", mode : "copy"

    input:
    tuple val(project), path(general_primer_flag), path(std_primer_flag)
    val outdir

    output:
    tuple val(project), path("*trimming_conductor_out.txt"), emit: trimming_conductor_out

    shell:

    '''
    gen_fwd_flag="$(sed '1q;d' !{general_primer_flag})"
    gen_rev_flag="$(sed '2q;d' !{general_primer_flag})"

    std_fwd_primer="$(grep 'F' !{std_primer_flag} || true)"
    std_rev_primer="$(grep 'R' !{std_primer_flag} || true)"

    fwd_trim_flag='none'
    rev_trim_flag='none'

    if [[ ! -z $std_fwd_primer ]]; then
        std_fwd_flag='1'
        fwd_trim_flag='std'
    else
        std_fwd_flag='0'
        if [[ $gen_fwd_flag -eq '1' ]]; then
            fwd_trim_flag='auto'
        else
            fwd_trim_flag='none'
        fi
    fi

    if [[ ! -z $std_rev_primer ]]; then
        std_rev_flag='1'
        rev_trim_flag='std'
    else
        std_rev_flag='0'
        if [[ $gen_rev_flag -eq '1' ]]; then
            rev_trim_flag='auto'
        else
            rev_trim_flag='none'
        fi
    fi

    echo "${fwd_trim_flag}\n${rev_trim_flag}" > trimming_conductor_out.txt
    '''

}