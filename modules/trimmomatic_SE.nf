
process TRIMMOMATIC_SE {
    
    label 'light'
    publishDir "${outdir}/${project}/${sampleId}/QC/", mode : "copy"

    input:
    tuple val(project), val(sampleId), path(fastq)
    val outdir

    output:
    tuple val(project), val(sampleId), path("*trimmed.fastq.gz"), emit: trimmed_fastq

    """
    trimmomatic SE -phred33 $fastq ${sampleId}.trimmed.fastq.gz LEADING:3 TRAILING:3 MINLEN:100 SLIDINGWINDOW:4:15
    """

}
