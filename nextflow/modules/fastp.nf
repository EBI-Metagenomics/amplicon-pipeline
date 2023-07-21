
process fastp {

    label 'light'
    // publishDir "${outdir}/merged/${project}", mode : "copy"
    container = '/hps/nobackup/rdf/metagenomics/singularity_cache/quay.io_microbiome-informatics_fastp:0.23.1.sif'

    input:
    tuple val(sampleId), path(fastq), val(project)
    val outdir

    output:
    tuple val(sampleId), path("*1_fastp.fastq.gz"), path("*2_fastp.fastq.gz"), val(project), emit: cleaned_fastq, optional: true

    """
    if [ -s ${fastq[0]} ] && [ -s ${fastq[1]} ]; then
        fastp -i ${fastq[0]} -I ${fastq[1]} --thread $task.cpus -o ${sampleId}_1_fastp.fastq.gz -O ${sampleId}_2_fastp.fastq.gz
    fi
    """

}
