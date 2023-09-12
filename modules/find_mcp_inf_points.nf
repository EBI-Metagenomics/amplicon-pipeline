
process FIND_MCP_INF_POINTS {
    // Find inflection points in conservation curves

    label 'light'

    input:
    tuple val(project), val(sampleId), val(var_region), path(mcp_cons_out)
    val outdir

    output:
    tuple val(project), val(sampleId), val(var_region), path("*inf_points.tsv"), emit: inf_points_out

    """
    if [[ -s ./$mcp_cons_out ]]; then
        python /hps/software/users/rdf/metagenomics/service-team/users/chrisata/asv_gen/bin/find_mcp_inflection_points_MERGED.py -i $mcp_cons_out -s res -o ./
    else
        touch res_inf_points.tsv
    fi
    """

}
