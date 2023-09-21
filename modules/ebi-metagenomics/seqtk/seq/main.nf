// This seqtk/seq module is simply copied from the already-existing nf-core module (https://nf-co.re/modules/seqtk_seq/, https://github.com/nf-core/modules/commit/726ee59cd9360a965d96ea9ea8770f16b8ddd6cc)
// This is because there are not currently any nf-core ways of adding modules from more than one nf-core repo

process SEQTK_SEQ {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::seqtk=1.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqtk:1.3--h5bf99c6_3' :
        'biocontainers/seqtk:1.3--h5bf99c6_3' }"

    input:
    tuple val(meta), path(fastx)

    output:
    tuple val(meta), path("*.gz")     , emit: fastx
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def extension = "fastq"
    if ("$fastx" ==~ /.+\.fasta|.+\.fasta.gz|.+\.fa|.+\.fa.gz|.+\.fas|.+\.fas.gz|.+\.fna|.+\.fna.gz/ || "$args" ==~ /\-[aA]/ ) {
        extension = "fasta"
    }
    """
    seqtk \\
        seq \\
        $args \\
        $fastx | \\
        gzip -c > ${prefix}.seqtk-seq.${extension}.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqtk: \$(echo \$(seqtk 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """
}
