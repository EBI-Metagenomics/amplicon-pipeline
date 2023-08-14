
process PARSE_CONDUCTOR {
    
    // Parse the trimming_conductor output and store 
    // the flags into environment variables 

    label 'light'
    // publishDir "${outdir}/merged/${project}", mode : "copy"

    input:
    tuple val(project), val(sampleId), path(trimming_conductor_out)
    val outdir

    output:
    tuple val(project), val(sampleId), env(fwd_flag), env(rev_flag), emit: conductor_out

    script:
    """
    fwd_flag=\$(sed '1q;d' "${trimming_conductor_out}")
    rev_flag=\$(sed '2q;d' "${trimming_conductor_out}")
    """
}
