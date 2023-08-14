
process ASSESS_MCP_CONS {
    // Use Most Common Prefix (MCP) method to generate curves of base conservation

    label 'light'

    input:
    val project
    val sampleId
    val fwd_flag
    val rev_flag
    path fastq
    val outdir

    output:
    tuple val(project), val(sampleId), path("*mcp_cons.tsv"), optional: true, emit: mcp_cons_out

    """
    if [[ ${fwd_flag} = "auto" ]] && [[ ${rev_flag} = "auto" ]]; then
        python /hps/software/users/rdf/metagenomics/service-team/users/chrisata/asv_gen/bin/assess_mcp_proportions_MERGED.py -i $fastq -s ${fastq.simpleName} -st FR -o ./
    elif [[ ${fwd_flag} = "auto" ]]; then
        python /hps/software/users/rdf/metagenomics/service-team/users/chrisata/asv_gen/bin/assess_mcp_proportions_MERGED.py -i $fastq -s ${fastq.simpleName} -st F -o ./
    elif [[ ${rev_flag} = "auto" ]]; then
        python /hps/software/users/rdf/metagenomics/service-team/users/chrisata/asv_gen/bin/assess_mcp_proportions_MERGED.py -i $fastq -s ${fastq.simpleName} -st R -o ./
    else
        touch ${fastq.simpleName}_mcp_cons.tsv
    fi
    """

}
