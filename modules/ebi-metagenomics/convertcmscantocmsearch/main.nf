
process CONVERTCMSCANTOCMSEARCH {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "microbiome-informatics/mgnify-pipelines-toolkit:1.0.2"

    input:
    tuple val(meta), path(cmscan_tblout)

    output:
    tuple val(meta), path("*cmsearch.tbl"), emit: cmsearch_tblout
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def is_compressed = cmscan_tblout.name.contains(".gz")
    def cmscan_input = is_compressed ? cmscan_tblout.name.replace(".gz", "") : cmscan_tblout
    """
    if [ "$is_compressed" == "true" ]; then
        gzip -c -d ${cmscan_tblout} > ${cmscan_input}
    fi
    convert_cmscan_to_cmsearch_tblout \\
        --input ${cmscan_input} \\
        --output ${meta.id}.cmsearch.tbl

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mgnify-pipelines-toolkit: \$(get_mpt_version)
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.cmsearch.tbl

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mgnify-pipelines-toolkit: \$(get_mpt_version)
    END_VERSIONS
    """
}
