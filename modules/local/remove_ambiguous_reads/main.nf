
process REMOVE_AMBIGUOUS_READS {
    // Run DADA2 pipeline including read-tracking

    label 'light'

    input:
    tuple val(meta), val(var_region), path(reads)

    output:
    tuple val(meta), val(var_region), path("*noambig*fastq.gz"), emit: noambig_out

    """
    if [[ ${meta.single_end} = true ]]; then
        python /hps/software/users/rdf/metagenomics/service-team/users/chrisata/asv_gen/bin/remove_ambiguous_reads.py -f $reads -s ${meta.id}
    else 
        python /hps/software/users/rdf/metagenomics/service-team/users/chrisata/asv_gen/bin/remove_ambiguous_reads.py -f ${reads[0]} -r ${reads[1]} -s ${meta.id}
    fi
    
    """

}
