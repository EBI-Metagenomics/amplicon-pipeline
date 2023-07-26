
process trimming_conductor {

    label 'light'
    publishDir "${outdir}/merged/${project}", mode : "copy"

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

    // case $comb_flag
    // in
    //     '1111')
    //         fwd_trim_flag='std'
    //         rev_trim_flag='std'
    //         ;;
        
    //     '0000')
    //         fwd_trim_flag='none'
    //         rev_trim_flag='none'
    //         ;;
    //     '1010')
    //         fwd_trim_flag='std'
    //         rev_trim_flag='none'
    //         ;;
    //     '0101')
    //         fwd_trim_flag='none'
    //         rev_trim_flag='std'
    //         ;;
    //     '1100')
    //         fwd_trim_flag='auto'
    //         rev_trim_flag='auto'
    //         ;;
    //     '1000')
    //         fwd_trim_flag='auto'
    //         rev_trim_flag='none'
    //         ;;
    //     '0100')
    //         fwd_trim_flag='none'
    //         rev_trim_flag='auto'
    //         ;;
    //     ''
        
    // esac


    // comb_flag=${gen_fwd_flag}${gen_rev_flag}${std_fwd_flag}${std_rev_flag}


    // if [ $comb_flag -eq '1111' ]; then
    //     trim_flag='true'
    // elif [ $comb_flag -eq '0000' ]; then
    //     trim_flag='false'
    // fi