
process FINAL_CONCAT_PRIMERS {

    label 'light'
    publishDir "${outdir}/${project}", mode : "copy"
    
    input:
    tuple val(project), val(sampleId), val(var_region), path(cat_primers)
    val outdir

    output:
    tuple val(project), val(sampleId), val(var_region), path("*final_concat_primers.fasta"), emit: final_concat_primers_out

    """
    cat *concat_primers.fasta > temp_concat_primers.fasta
    awk '!a[\$0]++' temp_concat_primers.fasta > ${sampleId}_final_concat_primers.fasta
    """

}
