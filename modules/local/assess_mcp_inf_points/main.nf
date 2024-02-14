
process ASSESS_MCP_INF_POINTS {
    // Select inflection points most likely to be primer cutoff points
    tag "$meta.id"
    label 'light'
    // publishDir "${outdir}/${project}/${sampleId}/primer-identification", mode : "copy" 

    input:
    tuple val(meta), path(inf_points_out), path(reads_merged)

    output:
    tuple val(meta), path("*auto_primers.fasta"), emit: auto_primer_out

    """
    if [[ -s ./$inf_points_out ]]; then
        python /hps/software/users/rdf/metagenomics/service-team/users/chrisata/asv_gen/bin/assess_inflection_point_mcp_MERGED.py -i $reads_merged -p $inf_points_out -s ${meta.id}_${meta.var_region} -o ./
    else
        touch ${meta.id}_${meta.var_region}_auto_primers.fasta
    fi
    """

}
