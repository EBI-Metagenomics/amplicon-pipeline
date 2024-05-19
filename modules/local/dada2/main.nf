
process DADA2 {
    // Run DADA2 pipeline including read-tracking
    tag "$meta.id"
    label "medium"
    // TODO: this env file doesn't exist
    conda "${moduleDir}/dada2_environment.yml"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*map.txt"), path("*asvs.fasta"), path("*_filt.fastq.gz"), emit: dada2_out
    path("*asv_counts.tsv"),                                                       emit: dada2_asv_counts
    tuple path("*chimeric.txt"), path("*matched.txt"),                             emit: dada2_stats
    
    script:
    if ( meta.single_end )
        """
        dada2.R ${meta.id} $reads
        """
    } else {
        """
        dada2.R ${meta.id} ${reads[0]} ${reads[1]}
        """
}
