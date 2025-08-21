
process EXTRACT_RRNA_HITS_FROM_UNMERGED_READS {
    tag "$meta.id"
    label 'process_medium'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqkit:2.9.0--h9ee0642_0' :
        'biocontainers/seqkit:2.9.0--h9ee0642_0' }"

    input:
    tuple val(meta), path(reads), path(cmsearch_deoverlap_out)

    output:
    tuple val(meta), path("*_extracted_read*.fastq.gz"), emit: extracted_reads
    path "versions.yml"                                , emit: versions

    script:
    """ 
    \$(cut -d' ' -f1 ${cmsearch_deoverlap_out} > ${meta.id}_extracted_ids.txt)

    if [[ ${meta.single_end} = true ]]; then
        seqkit grep -f ${meta.id}_extracted_ids.txt ${reads} -o ${meta.id}_extracted_reads.fastq.gz
    else
        seqkit grep -f ${meta.id}_extracted_ids.txt ${reads[0]} -o ${meta.id}_extracted_reads_1.fastq.gz
        seqkit grep -f ${meta.id}_extracted_ids.txt ${reads[1]} -o ${meta.id}_extracted_reads_2.fastq.gz
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$( seqkit | sed '3!d; s/Version: //' )
    END_VERSIONS
    """
}
