
process DADA2 {
    // Run DADA2 pipeline including read-tracking
    tag "$meta.id"
    label "dada2_resources"
    container 'docker://quay.io/microbiome-informatics/dada2:v1'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*map.txt"), path("*asvs.fasta"), path("*_filt.fastq.gz"), emit: dada2_out
    tuple val(meta), path("*_dada2_stats.tsv")                                     , emit: dada2_stats
    path "versions.yml"                                                            , emit: versions
    
    script:
    if ( meta.single_end ){
        """
        dada2.R ${meta.id} $reads

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            R: \$( R --version | head -1 | cut -d' ' -f3 )
            dada2: \$( R -e "suppressMessages(library(dada2));packageDescription('dada2')" | grep 'Version' | cut -d' ' -f2 )
        END_VERSIONS
        """
    } else {
        """
        dada2.R ${meta.id} ${reads[0]} ${reads[1]}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            R: \$( R --version | head -1 | cut -d' ' -f3 )
            dada2: \$( R -e "suppressMessages(library(dada2));packageDescription('dada2')" | grep 'Version' | cut -d' ' -f2 )
        END_VERSIONS
        """
    }

}
