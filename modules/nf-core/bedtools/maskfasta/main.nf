process BEDTOOLS_MASKFASTA {
    tag "$meta.id"
    label 'light'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.30.0--hc088bd4_0' :
        'biocontainers/bedtools:2.30.0--hc088bd4_0' }"

    input:
    tuple val(meta), path(fasta), path(bed)
    
    output:
    tuple val(meta), path("*.fa"), emit: fasta
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def fasta_unzipped = fasta.name.replace(".gz", "")
    """
    gzip -c -d $fasta > $fasta_unzipped
    bedtools \\
        maskfasta \\
        $args \\
        -fi $fasta_unzipped \\
        -bed $bed \\
        -fo ${prefix}.fa
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
    END_VERSIONS
    """
}
