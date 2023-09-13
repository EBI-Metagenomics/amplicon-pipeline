
process CLASSIFY_VAR_REGIONS {

    label 'light'
    publishDir "${outdir}/${project}/${sampleId}/amplified-region-inference", pattern : "*.concat.regions.txt" , mode : "copy"
    publishDir "${outdir}/${project}/${sampleId}/amplified-region-inference", pattern : "*.tsv" , mode : "copy"


    input:
    tuple val(project), val(sampleId), path(cmsearch_deoverlap_out)
    val outdir

    output:
    // tuple val(project), val(sampleId), path("*_regions.txt"), path("*seq_count.txt"), emit: classify_var_summary
    tuple val(project), val(sampleId), path ("*S.V*.txt"), optional: true, emit: classify_var_regions
    tuple val(project), val(sampleId), val("concat"), path("*.concat.regions.txt"), optional: true, emit: concat_var_regions
    tuple val(project), val(sampleId), path("*.tsv"), optional: true, emit: cdch_out

    """
    python /hps/software/users/rdf/metagenomics/service-team/users/chrisata/asv_gen/bin/classify_var_regions.py -d ./ -o ${cmsearch_deoverlap_out.simpleName} --statistics $cmsearch_deoverlap_out
    
    num_files="\$(ls *S.V*.txt | wc -l  )"
    if [ \$num_files -gt 1 ]; then
        cat *S.V*.txt > ${sampleId}.concat.regions.txt
    fi
    """

}
