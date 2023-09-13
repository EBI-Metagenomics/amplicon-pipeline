
process FASTP {
    // Run fastp on paired-end reads

    label 'light'
    container = '/hps/nobackup/rdf/metagenomics/singularity_cache/quay.io_microbiome-informatics_fastp:0.23.1.sif'

    input:
    tuple val(project), val(sampleId), path(fastq)
    val outdir

    output:
    tuple val(project), val(sampleId), path("*1_fastp.fastq.gz"), path("*2_fastp.fastq.gz"), emit: cleaned_fastq, optional: true

    """
    if [ -s ${fastq[0]} ] && [ -s ${fastq[1]} ]; then
        fastp -i ${fastq[0]} -I ${fastq[1]} --thread $task.cpus -o ${sampleId}_1_fastp.fastq.gz -O ${sampleId}_2_fastp.fastq.gz
    fi
    """

}
