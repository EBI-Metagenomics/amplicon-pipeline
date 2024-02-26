
process FIND_MCP_INF_POINTS {
    // Find inflection points in conservation curves
    tag "$meta.id"
    label 'light'
    conda "${projectDir}/conf/environment.yml"

    input:
    tuple val(meta), path(mcp_cons_out)

    output:
    tuple val(meta), path("*inf_points.tsv"), emit: inf_points_out

    """
    if [[ -s ./$mcp_cons_out ]]; then
        find_mcp_inflection_points -i $mcp_cons_out -s ${meta.id}_${meta.var_region} -o ./
    else
        touch ${meta.id}_${meta.var_region}_inf_points.tsv
    fi
    """

}
