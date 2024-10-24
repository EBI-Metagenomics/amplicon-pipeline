
process ASSESSMCPPROPORTIONS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "https://depot.galaxyproject.org/singularity/mgnify-pipelines-toolkit:${params.mpt_version}":
        "biocontainers/mgnify-pipelines-toolkit:${params.mpt_version}" }"

    input:
    tuple val(meta), val(fwd_flag), val(rev_flag), path(fastq)
    val library_check

    output:
    tuple val(meta), path("*.tsv") , emit: tsv
    tuple val(meta), env(check_out), optional: true, emit: library_check_out
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def var_region = "${meta.var_region}" ?: ""
    def assess_mcp_prop_prefix = "${prefix}_${var_region}"
    def strands = ""
    def library_check_input = "${assess_mcp_prop_prefix}_mcp_cons.tsv"

    if (fwd_flag == "auto" && rev_flag == "auto") {
        strands = "FR"
    } else if (fwd_flag == "auto") {
        strands = "F"
    } else if (rev_flag == "auto") {
        strands = "R"
    }

    if (strands == "") {
        """
        touch ${prefix}_${var_region}_mcp_cons.tsv

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            mgnify-pipelines-toolkit: \$(get_mpt_version)
        END_VERSIONS
        """
    } else if (library_check) {
        """
        assess_mcp_proportions \\
            -i ${fastq} \\
            -s ${assess_mcp_prop_prefix} \\
            -st ${strands} \\
            -o ./

        library_strategy_check \\
            -i ${library_check_input} \\
            -s ${prefix} \\
            -o ./

        check_out=\$(cat ${prefix}_library_check_out.txt)

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            mgnify-pipelines-toolkit: \$(get_mpt_version)
        END_VERSIONS
        """
    }  else {
        """
        assess_mcp_proportions \\
            -i ${fastq} \\
            -s ${assess_mcp_prop_prefix} \\
            -st ${strands} \\
            -o ./

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            mgnify-pipelines-toolkit: \$(get_mpt_version)
        END_VERSIONS
        """
    }

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def var_region = "${meta.var_region}"
    def assess_mcp_prop_prefix = "${prefix}_${var_region}"

    """
    touch ${assess_mcp_prop_prefix}_mcp_cons.tsv

    echo 'dummy' > ${assess_mcp_prop_prefix}_mcp_cons.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mgnify-pipelines-toolkit: \$(get_mpt_version)
    END_VERSIONS
    """
}
