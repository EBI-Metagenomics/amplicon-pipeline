
process EXTRACTCOORDS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "https://depot.galaxyproject.org/singularity/mgnify-pipelines-toolkit:${params.mpt_version}":
        "biocontainers/mgnify-pipelines-toolkit:${params.mpt_version}" }"

    input:
    tuple val(meta), path(easel_coords_fasta)
    tuple val(meta), path(matched_seqs_with_coords)

    output:
    tuple val(meta),  path("sequence-categorisation/*SSU.fasta")                         , optional: true, emit: ssu_fasta
    tuple val(meta),  path("sequence-categorisation/*LSU.fasta")                         , optional: true, emit: lsu_fasta
    tuple val(meta),  path("*concat_SSU_LSU_coords.txt")                                 , emit: concat_ssu_lsu_coords
    tuple val(meta),  path("sequence-categorisation/*rRNA_bacteria*.fasta")                 , optional: true, emit: rrna_bacteria
    tuple val(meta),  path("sequence-categorisation/*rRNA_archaea*.fasta")                  , optional: true, emit: rrna_archaea
    tuple val(meta),  path("sequence-categorisation/*rRNA_eukarya*.fasta")                  , optional: true, emit: eukarya
    path "versions.yml"                                                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    get_subunits -i $easel_coords_fasta -n ${prefix} --separate-subunits-by-models
    get_subunits_coords -i $matched_seqs_with_coords -s SSU -l LSU
    cat SSU_coords LSU_coords > ${prefix}_concat_SSU_LSU_coords.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        extractcoords: 0.1.2
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}_concat_SSU_LSU_coords.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        extractcoords: 0.1.2
    END_VERSIONS
    """
}
