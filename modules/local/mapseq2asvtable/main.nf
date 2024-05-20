
process MAPSEQ2ASVTABLE {
    tag "$meta.id"
    label 'very_light' // Will likely need to give this task more CPUs 

    // TODO: add a container please

    input:
    tuple val(meta), path(mapseq_out)
    val(db_label)

    output:
    tuple val(meta), path("*.tsv"), emit: asvtaxtable

    """
    mapseq2asvtable.py -i $mapseq_out -l $db_label -s ${meta.id}
    """

}
