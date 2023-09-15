
process MAPSEQ2BIOM {

    label 'light' // Will likely need to give this task more CPUs
    publishDir "${outdir}/${project}/${sampleId}/taxonomy-summary/${label}", mode : "copy"
 
    
    input:
    tuple val(project), val(sampleId), path(mapseq_out)
    tuple path(db_fasta), path(db_tax), path(db_otu), path(db_mscluster), val(label)
    val outdir

    output:
    tuple val(project), val(sampleId), path("${mapseq_out.simpleName}.tsv"), path("${mapseq_out.simpleName}.txt"), path("*notaxid.tsv"), emit: mapseq2biom_out

    """
    python /hps/software/users/rdf/metagenomics/service-team/users/chrisata/asv_gen/bin/mapseq2biom.py \
        --out-file ${mapseq_out.simpleName}.tsv \
        --krona ${mapseq_out.simpleName}.txt \
        --no-tax-id-file ${mapseq_out.simpleName}.notaxid.tsv \
        --taxid \
        --label ${label} \
        --query ${mapseq_out} \
        --otu-table ${db_otu}
    """

}
