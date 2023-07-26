
process parse_conductor {

    label 'light'
    // publishDir "${outdir}/merged/${project}", mode : "copy"

    input:
    tuple val(project), path(trimming_conductor_out)
    val outdir

    output:
    tuple val(project), env(fwd_flag), env(rev_flag), emit: trimming_flags_out

    script:
    """
    fwd_flag=\$(sed '1q;d' "${trimming_conductor_out}")
    rev_flag=\$(sed '2q;d' "${trimming_conductor_out}")
    """
}