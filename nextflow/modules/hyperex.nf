
process hyperex {

    publishDir "${outdir}/merged/${project}", mode : "copy"

    input:
    tuple path(fasta), val(project)
    val outdir

    output: // TODO: Probably don't need hyperex.gff, keep it for now
    tuple path("*_hyperex.fa"), path("*_hyperex.gff"), val(project), emit: hyperex_out

    """
    /hps/nobackup/rdf/metagenomics/service-team/users/chrisata/hyperex/hyperex-0.1.0/hyperex -p ${fasta.simpleName}_hyperex $fasta
    """
}
