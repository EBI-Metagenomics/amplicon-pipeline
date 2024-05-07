
process CMSEARCHTBLOUTDEOVERLAP {
    tag "$meta.id"
    label 'very_light'

    conda "bioconda::cmsearch_tblout_deoverlap=0.09"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cmsearch_tblout_deoverlap:0.09--pl5321hdfd78af_0':
        'biocontainers/cmsearch_tblout_deoverlap:0.09--pl5321hdfd78af_0' }"

    input:
    tuple val(meta), path(cmsearch_tblout)
    path(clanin)

    output:
    tuple val(meta), path("*.deoverlapped"), emit: cmsearch_tblout_deoverlapped
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def cmsearch_tblout_unzipped = cmsearch_tblout.name.replace(".gz", "")
    """
    gzip -c -d $cmsearch_tblout > $cmsearch_tblout_unzipped
    cmsearch-deoverlap.pl $args \
        --clanin $clanin \
        $cmsearch_tblout_unzipped

    mv ${cmsearch_tblout_unzipped}.deoverlapped ${prefix}.tblout.deoverlapped

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cmsearch_tblout_deoverlap: \$(echo \$(cmsearch-deoverlap.pl 2>&1) | grep -o 'cmsearch-deoverlap v[0-9]\\+\\.[0-9]\\+' | sed 's/cmsearch-deoverlap //g' )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tblout.deoverlapped

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cmsearch_tblout_deoverlap: \$(echo \$(cmsearch-deoverlap.pl 2>&1) | grep -o 'cmsearch-deoverlap v[0-9]\\+\\.[0-9]\\+' | sed 's/cmsearch-deoverlap //g' )
    END_VERSIONS
    """
}
