
process std_primer_flag {

    label 'light'
    publishDir "${outdir}/merged/${project}", mode : "copy"

    input:
    tuple path(fastq), val(project)
    val outdir

    output:
    path "*std_primer_out.txt", emit: std_primer_out

    """
    python /hps/software/users/rdf/metagenomics/service-team/users/chrisata/asv_gen/bin/standard_primer_agrep.py -i $fastq -s ${fastq.simpleName} -o ./
    """

}
