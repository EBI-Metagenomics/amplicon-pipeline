
process INFERNAL_CMSEARCH {

    tag "$meta.id"

    label 'process_low'

    conda "bioconda::infernal=1.1.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/infernal:1.1.4--pl5321hec16e2b_1':
        'biocontainers/infernal:1.1.4--pl5321hec16e2b_1' }"

    input:
    tuple val(meta), path(seqdb)
    path covariance_model_database

    output:
    tuple val(meta), path("*.cmsearch_matches.tbl.gz"), emit: cmsearch_tbl
    path "versions.yml"                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def is_compressed = seqdb.name.endsWith(".gz")
    def seqdb_name = seqdb.name.replace(".gz", "")
    """
    if [ "$is_compressed" == "true" ]; then
        gzip -c -d $seqdb > $seqdb_name
    fi

    cmsearch \\
        --cpu $task.cpus \\
        --noali \\
        --hmmonly \\
        $args \\
        -Z 1000 \\
        -o /dev/null \\
        --tblout ${prefix}.cmsearch_matches.tbl \\
        $covariance_model_database \\
        $seqdb_name

    gzip -n ${prefix}.cmsearch_matches.tbl

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cmsearch: \$(cmsearch -h | grep -o '^# INFERNAL [0-9.]*' | sed 's/^# INFERNAL *//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.cmsearch_matches.tbl
    gzip ${prefix}.cmsearch_matches.tbl

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        infernal: \$(cmsearch -h | grep -o '^# INFERNAL [0-9.]*' | sed 's/^# INFERNAL *//')
    END_VERSIONS
    """
}
