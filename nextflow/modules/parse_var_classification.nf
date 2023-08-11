
process PARSE_VAR_CLASSIFICATION {

    label 'light'

    input:
    tuple val(project), path(classify_var_out), path(classify_var_seqcount)
    path amplified_regions
    val outdir

    output:
    stdout
    // tuple val(project), path("*.tsv"), path("*seq_count.txt"), emit: classify_var_summary
    // path "*.V*.txt", optional: true, emit: classify_var_regions

    """
    echo ${amplified_regions[0]}
    echo ${amplified_regions[1]}
    cat $classify_var_out 

    """
    // amp_regions=\$(cut -f4-5 $classify_var_out | tail -n +2 | sed 's/\\t/;/g')
    // echo \$amp_regions
}
