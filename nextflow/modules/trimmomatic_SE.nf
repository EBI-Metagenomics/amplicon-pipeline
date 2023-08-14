
process TRIMMOMATIC_SE {
    
    label 'light'

    input:
    tuple val(project), val(sampleId), path(fastq)
    val outdir

    output:
    tuple val(project), val(sampleId), path("*trimmed.fastq.gz"), emit: trimmed_fastq

    """
    trimmomatic SE -phred33 $fastq ${fastq.simpleName}.trimmed.fastq.gz LEADING:3 TRAILING:3 MINLEN:100 SLIDINGWINDOW:4:15
    """

}
