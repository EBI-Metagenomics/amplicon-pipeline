
process MAPSEQ2ASVTABLE {
    tag "$meta.id"
    label 'very_light' // Will likely need to give this task more CPUs 
    conda "${projectDir}/conf/environment.yml"
    // TODO: uncomment container when you release fix to mpt
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     "https://depot.galaxyproject.org/singularity/mgnify-pipelines-toolkit:${params.mpt_version}":
    //     "biocontainers/mgnify-pipelines-toolkit:${params.mpt_version}" }"

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
