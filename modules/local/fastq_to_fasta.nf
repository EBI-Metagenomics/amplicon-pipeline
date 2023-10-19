
process FASTQ_TO_FASTA {
    // Convert fastq to fasta
    
    label 'light'
    publishDir "${outdir}/${project}/${sampleId}/QC/", mode : "copy"
    container = '/hps/nobackup/rdf/metagenomics/service-team/singularity-cache/quay.io_biocontainers_seqtk:1.3.sif'

    input:
    tuple val(project), val(sampleId), path(fastq)
    val outdir

    output:
    tuple val(project), val(sampleId), path("*.fasta"), emit: merged_fasta

    """
    seqtk seq -a $fastq > ${sampleId}.fasta
    """

}
