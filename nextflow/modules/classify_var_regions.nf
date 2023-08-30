
process CLASSIFY_VAR_REGIONS {

    label 'light'

    input:
    tuple val(project), val(sampleId), path(cmsearch_deoverlap_out)
    val outdir

    output:
    tuple val(project), val(sampleId), path("*_regions.txt"), path("*seq_count.txt"), emit: classify_var_summary
    path "*S.V*.txt", optional: true, emit: classify_var_regions

    """
    echo "meow"
    python /hps/software/users/rdf/metagenomics/service-team/users/chrisata/asv_gen/bin/classify_var_regions.py -d ./ -o ${cmsearch_deoverlap_out.simpleName} --statistics $cmsearch_deoverlap_out
    """

}
