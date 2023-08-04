
process std_primer_flag {
    // Check for presence of standard library of primers (stored in ./data/standard_primers)

    label 'light'
    publishDir "${outdir}/merged/${project}", mode : "copy"

    input:
    tuple  val(project), path(fastq)
    val outdir

    output:
    tuple val(project), path("*std_primers.fasta"), emit: std_primer_out

    """
    python /hps/software/users/rdf/metagenomics/service-team/users/chrisata/asv_gen/bin/standard_primer_agrep.py -i $fastq -s ${fastq.simpleName} -o ./
    """

}
