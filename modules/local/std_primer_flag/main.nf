
process STD_PRIMER_FLAG {
    // Check for presence of standard library of primers (stored in ./data/standard_primers)
    tag "$meta.id"
    label 'very_light'
    container "oras://community.wave.seqera.io/library/mi-pimento:0.0.4--1795182b3e13ee89"

    input:
    tuple val(meta), path(reads_merged)
    path(std_primer_library)

    output:
    tuple val(meta), path("*std_primers.fasta"), emit: std_primer_out
    path "*std_primer_out.txt"
    path "versions.yml"                        , emit: versions


    script:
    """
    pimento std -i ${reads_merged} -p ${std_primer_library} -o ${meta.id}_${meta.var_region}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mi-pimento: \$( pimento --version | cut -d" " -f3 )
    END_VERSIONS
    """

}
