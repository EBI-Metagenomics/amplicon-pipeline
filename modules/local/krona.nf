
process KRONA {
    tag "$meta.id"
    label 'light'
    // publishDir "${outdir}/${project}/${sampleId}/taxonomy-summary/${label}", mode : "copy"
    container = '/hps/nobackup/rdf/metagenomics/service-team/singularity-cache/quay.io-biocontainers-krona-2.7.1--pl5321hdfd78af_7.img'

    input:
    tuple val(meta), path(otu_counts)
    tuple path(db_fasta), path(db_tax), path(db_otu), path(db_mscluster), val(label)

    output:
    tuple val(meta), path('*krona.html'), emit: krona_out
    
    script:
    prefix = otu_counts.baseName

    """
    ktImportText -o ${prefix}_krona.html $otu_counts
    """

}
