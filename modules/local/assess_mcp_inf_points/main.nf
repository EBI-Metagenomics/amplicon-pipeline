
process ASSESS_MCP_INF_POINTS {
    // Select inflection points most likely to be primer cutoff points
    tag "$meta.id"
    label 'light'
    conda "bioconda::mgnify-pipelines-toolkit=0.1.0"

    input:
    tuple val(meta), path(inf_points_out), path(reads_merged)

    output:
    tuple val(meta), path("*auto_primers.fasta"), emit: auto_primer_out

    """
    if [[ -s ./$inf_points_out ]]; then
        assess_inflection_point_mcp -i $reads_merged -p $inf_points_out -s ${meta.id}_${meta.var_region} -o ./
    else
        touch ${meta.id}_${meta.var_region}_auto_primers.fasta
    fi
    """

}
