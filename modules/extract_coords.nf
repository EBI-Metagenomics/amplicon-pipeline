
process EXTRACT_COORDS {

    label 'light'

    input:
    tuple val(project), val(sampleId), path(fasta)
    path(matched_seqs_with_coords)
    val outdir

    output:
    tuple val(project), val(sampleId)
    tuple val(project), val(sampleId), path("sequence-categorisation/*SSU.fasta"), optional: true, emit: ssu_fasta
    tuple val(project), val(sampleId), path("sequence-categorisation/*LSU.fasta"), optional: true, emit: lsu_fasta
    tuple val(project), val(sampleId), path("*concat_SSU_LSU_coords"), optional: true, emit: concat_ssu_lsu_coords
    
    """
    python /hps/software/users/rdf/metagenomics/service-team/users/chrisata/asv_gen/bin/get_subunits.py -i $fasta -n $sampleId
    python /hps/software/users/rdf/metagenomics/service-team/users/chrisata/asv_gen/bin/get_subunits_coords.py -i $matched_seqs_with_coords -s SSU -l LSU
    cat SSU_coords LSU_coords > ${sampleId}_concat_SSU_LSU_coords
    """

}
