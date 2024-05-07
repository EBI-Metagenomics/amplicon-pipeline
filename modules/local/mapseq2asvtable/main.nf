
process MAPSEQ2ASVTABLE {
    tag "$meta.id"
    // conda "${projectDir}/conf/environment.yml"
    label 'very_light' // Will likely need to give this task more CPUs 
    
    input:
    tuple val(meta), path(mapseq_out)
    val(db_label)

    output:
    tuple val(meta), path("*.tsv"), emit: asvtaxtable

    """
    python /hps/software/users/rdf/metagenomics/service-team/users/chrisata/asv_gen/bin/mapseq2asvtable.py -i $mapseq_out -l $db_label -s ${meta.id}
    """

}
