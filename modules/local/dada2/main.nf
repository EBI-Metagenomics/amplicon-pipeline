
process DADA2 {
    // Run DADA2 pipeline including read-tracking
    tag "$meta.id"
    label "medium"

    input:
    tuple val(meta), path(reads)
    path dada2_db
    val(db_label)

    output:
    tuple val(meta), path("*map.txt"), path("*asvs.fasta"), path("*_filt.fastq.gz"), optional: true, emit: dada2_out
    path("*asv_counts.tsv"), emit: dada2_asv_counts
    tuple path("*chimeric.txt"), path("*matched.txt"), emit: dada2_stats
    
    """
    echo "meoww"
    if [[ ${meta.single_end} = true ]]; then
        Rscript /hps/software/users/rdf/metagenomics/service-team/users/chrisata/asv_gen/bin/dada2.R ${meta.id} $db_label $dada2_db $reads
    else
        Rscript /hps/software/users/rdf/metagenomics/service-team/users/chrisata/asv_gen/bin/dada2.R ${meta.id} $db_label $dada2_db ${reads[0]} ${reads[1]}
    fi
    """

}
