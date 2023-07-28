
process assess_mcp_inf_points {

    label 'light'
    publishDir "${outdir}/merged/${project}", mode : "copy"

    input:
    tuple val(project), path(inf_points_out), path(fastq)
    val outdir

    output:
    tuple val(project), path("*auto_primers.fasta"), emit: auto_primer_out

    """
    python /hps/software/users/rdf/metagenomics/service-team/users/chrisata/asv_gen/bin/assess_inflection_point_mcp_MERGED.py -i $fastq -p $inf_points_out -s ${fastq.simpleName} -o ./
    """

}
