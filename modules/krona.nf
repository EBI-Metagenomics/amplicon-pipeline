
process KRONA {

    label 'light'
    publishDir "${outdir}/${project}/${sampleId}/taxonomy-summary/${label}", mode : "copy"
    container = '/hps/nobackup/rdf/metagenomics/singularity_cache_nextflow/quay.io-biocontainers-krona-2.7.1--pl5321hdfd78af_7.img'

    input:
    tuple val(project), val(sampleId), path(otu_counts)
    tuple path(db_fasta), path(db_tax), path(db_otu), path(db_mscluster), val(label)
    val outdir

    output:
    tuple val(project), val(sampleId), path('*krona.html'), emit: krona_out

    """
    ktImportText -o ${otu_counts.baseName}_krona.html $otu_counts
    """

}
