
process GENERAL_PRIMER_FLAG {
    // Check for the presence of primers in general
    tag "$meta.id"
    label 'very_light'
    container "oras://community.wave.seqera.io/library/mi-pimento:0.0.4--1795182b3e13ee89"

    input:
    tuple val(meta), path(reads_merged)

    output:
    tuple val(meta), path("*general_primer_out.txt"), emit: general_primer_out
    path "versions.yml"                             , emit: versions

    script:
    """
    pimento are_there_primers -i ${reads_merged} -o ${meta.id}_${meta.var_region}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mi-pimento: \$( pimento --version | cut -d" " -f3 )
    END_VERSIONS
    """
}
