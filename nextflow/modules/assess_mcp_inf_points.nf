
process assess_mcp_inf_points {

    label 'light'
    publishDir "${outdir}/merged/${project}", mode : "copy"

    input:
    tuple val(project), path(inf_points_out), path(fastq)
    val outdir

    output:
    path "*cutoff.txt"

    """
    python /hps/software/users/rdf/metagenomics/service-team/users/chrisata/asv_gen/bin/assess_inflection_point_mcp_MERGED.py -i $fastq -p $inf_points_out -s ${fastq.simpleName} -o ./
    """

}
