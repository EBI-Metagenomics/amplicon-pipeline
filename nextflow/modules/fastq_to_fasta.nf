
process fastq_to_fasta {
    
    label 'light'
    // publishDir "${outdir}/merged/${project}", mode : "copy"
    container = '/hps/nobackup/rdf/metagenomics/singularity_cache/quay.io_biocontainers_seqtk:1.3.sif'

    input:
    tuple path(fastq), val(project)
    val outdir

    output:
    tuple path("*.fasta"), val(project), emit: conv_fasta

    """
    seqtk seq -a $fastq > ${fastq.simpleName}.fasta
    """

}
