
process REMOVE_AMBIGUOUS_READS {
    // Run DADA2 pipeline including read-tracking
    tag "$meta.id"
    label 'process_low'
    conda "${projectDir}/conf/environment.yml"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*noambig*fastq.gz"), emit: noambig_out

    """
    if [[ ${meta.single_end} = true ]]; then
        remove_ambiguous_reads -f $reads -s ${meta.id}
    else 
        remove_ambiguous_reads -f ${reads[0]} -r ${reads[1]} -s ${meta.id}
    fi
    """

}
