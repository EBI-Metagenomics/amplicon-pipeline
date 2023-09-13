
process ASSESS_MCP_INF_POINTS {
    // Select inflection points most likely to be primer cutoff points

    label 'light'
    publishDir "${outdir}/${project}/${sampleId}/primer-identification", mode : "copy" 

    input:
    tuple val(project), val(sampleId), val(var_region), path(inf_points_out), path(fastq)
    val outdir

    output:
    tuple val(project), val(sampleId), val(var_region), path("*auto_primers.fasta"), emit: auto_primer_out

    """
    if [[ -s ./$inf_points_out ]]; then
        python /hps/software/users/rdf/metagenomics/service-team/users/chrisata/asv_gen/bin/assess_inflection_point_mcp_MERGED.py -i $fastq -p $inf_points_out -s ${fastq.simpleName}_${var_region} -o ./
    else
        touch ${fastq.simpleName}_${var_region}_auto_primers.fasta
    fi
    """

}
