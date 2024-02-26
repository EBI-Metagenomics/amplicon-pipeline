
process MAPSEQ2BIOM {
    tag "$meta.id"
    conda "${projectDir}/conf/environment.yml"
    label 'light' // Will likely need to give this task more CPUs 
    
    input:
    tuple val(meta), path(mapseq_out)
    tuple path(db_fasta), path(db_tax), path(db_otu), path(db_mscluster), val(label)

    output:
    tuple val(meta), path("${meta.id}.txt"), emit: mapseq2biom_krona_out
    tuple val(meta), path("${meta.id}.tsv"), path("*notaxid.tsv"), emit: mapseq2biom_out

    """
    mapseq2biom \
        --out-file ${meta.id}.tsv \
        --krona ${meta.id}.txt \
        --no-tax-id-file ${meta.id}.notaxid.tsv \
        --taxid \
        --label ${label} \
        --query ${mapseq_out} \
        --otu-table ${db_otu}
    """

}
