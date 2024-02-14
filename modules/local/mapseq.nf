
process MAPSEQ {

    label 'process_medium' // Will likely need to give this task more CPUs 
    container = '/hps/nobackup/rdf/metagenomics/service-team/singularity-cache/quay.io-biocontainers-mapseq-2.1.1--ha34dc8c_0.img'

    input:
    tuple val(meta), path(subunit_fasta)
    tuple path(db_fasta), path(db_tax), path(db_otu), path(db_mscluster), val(label)

    output:
    tuple val(meta), path('*.mseq'), emit: mapseq_out

    """
    mapseq $subunit_fasta $db_fasta $db_tax -tophits 80 -topotus 40 -outfmt 'simple' > ${meta.id}.mseq
    """

}
