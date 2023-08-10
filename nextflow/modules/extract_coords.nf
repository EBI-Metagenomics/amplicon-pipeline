
process extract_coords {

    label 'light'

    input:
    tuple val(project), path(fasta)
    val outdir

    output:
    val project
    path "sequence-categorisation/*SSU.fasta", optional: true, emit: ssu_fasta
    path "sequence-categorisation/*LSU.fasta", optional: true, emit: lsu_fasta
    
    """
    python /hps/software/users/rdf/metagenomics/service-team/users/chrisata/asv_gen/bin/get_subunits.py -i $fasta -n ${fasta.simpleName}
    """

}
