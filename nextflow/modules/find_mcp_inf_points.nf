
process find_mcp_inf_points {
    // Find inflection points in conservation curves

    label 'light'

    input:
    tuple val(project), path(mcp_cons_out)
    val outdir

    output:
    tuple val(project), path("*inf_points.tsv"), emit: inf_points_out

    """
    if [[ -s ./$mcp_cons_out ]]; then
        python /hps/software/users/rdf/metagenomics/service-team/users/chrisata/asv_gen/bin/find_mcp_inflection_points_MERGED.py -i $mcp_cons_out -s res -o ./
    else
        touch res_inf_points.tsv
    fi
    """

}
