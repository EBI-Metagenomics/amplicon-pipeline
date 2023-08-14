
process EXTRACT_COORDS {

    label 'light'

    input:
    tuple val(project), val(sampleId), path(fasta)
    val outdir

    output:
    tuple val(project), val(sampleId)
    tuple val(project), val(sampleId), path("sequence-categorisation/*SSU.fasta"), optional: true, emit: ssu_fasta
    tuple val(project), val(sampleId), path("sequence-categorisation/*LSU.fasta"), optional: true, emit: lsu_fasta
    
    """
    python /hps/software/users/rdf/metagenomics/service-team/users/chrisata/asv_gen/bin/get_subunits.py -i $fasta -n ${fasta.simpleName}
    """

}
