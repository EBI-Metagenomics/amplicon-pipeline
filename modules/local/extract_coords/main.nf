
process EXTRACT_COORDS {
    tag "$meta.id"
    label 'light'
    conda "${projectDir}/conf/environment.yml"

    input:
    tuple val(meta), path(easel_coords)
    tuple val(meta), path(matched_seqs_with_coords)

    output:
    tuple val(meta),  path("sequence-categorisation/*SSU.fasta"), optional: true, emit: ssu_fasta
    tuple val(meta),  path("sequence-categorisation/*LSU.fasta"), optional: true, emit: lsu_fasta
    tuple val(meta),  path("*concat_SSU_LSU_coords"), optional: true, emit: concat_ssu_lsu_coords


    """
    get_subunits -i $easel_coords -n ${meta.id}
    get_subunits_coords -i $matched_seqs_with_coords -s SSU -l LSU
    cat SSU_coords LSU_coords > ${meta.id}_concat_SSU_LSU_coords
    """

}
