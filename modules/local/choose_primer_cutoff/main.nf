
process CHOOSE_PRIMER_CUTOFF {
    // Select inflection points most likely to be primer cutoff points
    tag "$meta.id"
    label 'very_light'
    container "oras://community.wave.seqera.io/library/mi-pimento:0.0.4--1795182b3e13ee89"

    input:
    tuple val(meta), path(cutoffs), path(reads_merged)

    output:
    tuple val(meta), path("*auto_primers.fasta"), emit: auto_primer_out
    path "versions.yml"                         , emit: versions

    script:
    """
    if [[ -s ./$cutoffs ]]; then
        pimento choose_primer_cutoff -i ${reads_merged} -p ${cutoffs} -o ${meta.id}_${meta.var_region}
    else
        touch ${meta.id}_${meta.var_region}_auto_primers.fasta
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mi-pimento: \$( pimento --version | cut -d" " -f3 )
    END_VERSIONS
    """

}
