process SEQTK_SUBSEQ {
    tag "$sequences"
    label 'very_light'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqtk:1.3--h5bf99c6_3' :
        'biocontainers/seqtk:1.3--h5bf99c6_3' }"

    input:
    tuple val(meta), path(filter_list), path(sequences)

    output:
    tuple val(meta), val(var_region), path("*.gz"), emit: extracted_var_out
    tuple val(meta), val(var_region), path(filter_list), emit: extracted_var_path
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def ext = "fasta"
    var_region = "${filter_list.baseName.split('\\.')[1,2].join('-')}"
    if ("$sequences" ==~ /.+\.fq|.+\.fq.gz|.+\.fastq|.+\.fastq.gz/) {
        ext = "fastq"
    }
    
    """
    seqtk \\
        subseq \\
        $args \\
        $sequences \\
        $filter_list | \\
        gzip --no-name > ${prefix}_${var_region}_extracted.${ext}.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqtk: \$(echo \$(seqtk 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """
}
