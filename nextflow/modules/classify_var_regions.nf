
process classify_var_regions {

    label 'light'

    input:
    tuple val(project), path(cmsearch_deoverlap_out)
    val outdir

    output:
    tuple val(project), path("*_regions.txt"), path("*seq_count.txt"), emit: classify_var_summary
    path "*.V*.txt", optional: true, emit: classify_var_regions

    """
    python /hps/software/users/rdf/metagenomics/service-team/users/chrisata/asv_gen/bin/classify_var_regions.py -d ./ -o ${cmsearch_deoverlap_out.simpleName} --statistics $cmsearch_deoverlap_out
    """

}
