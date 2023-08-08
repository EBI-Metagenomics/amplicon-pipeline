
process assess_mcp_inf_points {
    // Select inflection points most likely to be primer cutoff points

    label 'light'
    publishDir "${outdir}/${project}", mode : "copy"

    input:
    tuple val(project), path(inf_points_out), path(fastq)
    val outdir

    output:
    tuple val(project), path("*auto_primers.fasta"), emit: auto_primer_out

    """
    if [[ -s ./$inf_points_out ]]; then
        python /hps/software/users/rdf/metagenomics/service-team/users/chrisata/asv_gen/bin/assess_inflection_point_mcp_MERGED.py -i $fastq -p $inf_points_out -s ${fastq.simpleName} -o ./
    else
        touch ${fastq.simpleName}_auto_primers.fasta
    fi
    """

}
