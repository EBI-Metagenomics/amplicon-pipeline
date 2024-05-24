
process STD_PRIMER_FLAG {
    // Check for presence of standard library of primers (stored in ./data/standard_primers)
    tag "$meta.id"
    label 'very_light'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "https://depot.galaxyproject.org/singularity/mgnify-pipelines-toolkit:${params.mpt_version}":
        "biocontainers/mgnify-pipelines-toolkit:${params.mpt_version}" }"

    input:
    tuple val(meta), path(reads_merged)
    path(std_primer_library)

    output:
    tuple val(meta), path("*std_primers.fasta"), emit: std_primer_out
    path "*std_primer_out.txt"
    path "versions.yml"                        , emit: versions


    script:
    """
    standard_primer_matching -i $reads_merged -p $std_primer_library -s ${meta.id}_${meta.var_region} -o ./

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mgnify-pipelines-toolkit: ${params.mpt_version}
    END_VERSIONS
    """

}
