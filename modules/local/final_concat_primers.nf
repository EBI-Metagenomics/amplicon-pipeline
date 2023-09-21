
process FINAL_CONCAT_PRIMERS {

    label 'light'
    // publishDir "${outdir}/${project}/${sampleId}/primer-identification", mode : "copy" 
    
    input:
    tuple val(meta), val(var_region), path(cat_primers)

    output:
    tuple val(meta), val(var_region), path("*final_concat_primers.fasta"), emit: final_concat_primers_out

    """
    cat *concat_primers.fasta > temp_concat_primers.fasta
    awk '!a[\$0]++' temp_concat_primers.fasta > ${meta.id}_final_concat_primers.fasta
    """

}
