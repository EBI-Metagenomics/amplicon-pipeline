
process REMOVE_AMBIGUOUS_READS {
    // Run DADA2 pipeline including read-tracking

    label 'light'

    input:
    tuple val(project), val(sampleId), val(var_region), path(fastq_1), path(fastq_2)
    val outdir

    output:
    tuple val(project), val(sampleId), val(var_region), path("*noambig_1.fastq.gz"), path("*noambig_2.fastq.gz"), emit: noambig_out

    """
    python /hps/software/users/rdf/metagenomics/service-team/users/chrisata/asv_gen/bin/remove_ambiguous_reads.py -f $fastq_1 -r $fastq_2 -s $sampleId
    """

}
