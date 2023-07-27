
process find_mcp_inf_points {

    label 'light'
    // publishDir "${outdir}/merged/${project}", mode : "copy"

    input:
    tuple val(project), path(mcp_cons_out)
    val outdir

    output:
    tuple val(project), path("*inf_points.tsv"), emit: inf_points_out

    """
    python /hps/software/users/rdf/metagenomics/service-team/users/chrisata/asv_gen/bin/find_mcp_inflection_points_MERGED.py -i $mcp_cons_out -s res -o ./

    """

}
