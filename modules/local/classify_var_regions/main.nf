
process CLASSIFY_VAR_REGIONS {
    tag "$meta.id"
    label 'light'
    // publishDir "${outdir}/${project}/${sampleId}/amplified-region-inference", pattern : "*.concat.regions.txt" , mode : "copy"
    // publishDir "${outdir}/${project}/${sampleId}/amplified-region-inference", pattern : "*.tsv" , mode : "copy"


    input:
    tuple val(meta), path(cmsearch_deoverlap_out)

    output:
    tuple val(meta), path ("*S.V*.txt"), optional: true, emit: classify_var_regions
    tuple val(meta), val("concat"), path("*.concat.regions.txt"), optional: true, emit: concat_var_regions
    tuple val(meta), path("*.tsv"), optional: true, emit: cdch_out

    """
    python /hps/software/users/rdf/metagenomics/service-team/users/chrisata/asv_gen/bin/classify_var_regions.py -d ./ -o ${meta.id} --statistics $cmsearch_deoverlap_out
    
    num_files="\$(ls *S.V*.txt | wc -l  )"
    if [ \$num_files -gt 1 ]; then
        cat *S.V*.txt > ${meta.id}.concat.regions.txt
    fi
    
    """

}
