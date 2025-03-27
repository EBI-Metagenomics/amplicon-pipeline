
process FIND_PRIMER_CUTOFFS {
    // Find inflection points in conservation curves
    tag "$meta.id"
    label 'very_light'
    container "oras://community.wave.seqera.io/library/mi-pimento:0.0.4--1795182b3e13ee89"

    input:
    tuple val(meta), path(bcv)

    output:
    tuple val(meta), path("*cutoffs.tsv"), emit: cutoffs
    path "versions.yml"                  , emit: versions

    script:
    """
    if [[ -s ./$bcv ]]; then
        pimento find_cutoffs -i ${bcv} -o ${meta.id}_${meta.var_region}
    else
        touch ${meta.id}_${meta.var_region}_cutoffs.tsv
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mi-pimento: \$( pimento --version | cut -d" " -f3 )
    END_VERSIONS
    """
}
