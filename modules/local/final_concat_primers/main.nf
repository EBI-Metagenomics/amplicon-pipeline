
process FINAL_CONCAT_PRIMERS {
    tag "$meta.id"
    label 'light'
    
    input:
    tuple val(meta), path(cat_primers)

    output:
    tuple val(meta), path("*final_concat_primers.fasta"), optional:true, emit: final_concat_primers_out

    """
    cat *concat_primers.fasta > temp_concat_primers.fasta
    if [[ -s temp_concat_primers.fasta ]]; then
        awk '!a[\$0]++' temp_concat_primers.fasta > ${meta.id}_final_concat_primers.fasta
    else
       touch ${meta.id}_final_concat_primers.fasta
    fi
    """

}
