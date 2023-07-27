
process seqprep_merge {

    label 'light'
    publishDir "${outdir}/merged/${project}", mode : "copy"
    container = '/hps/nobackup/rdf/metagenomics/singularity_cache_nextflow/quay.io-biocontainers-seqprep-1.3.2--hed695b0_4.img'

    input:
    tuple val(sampleId), path(fastq_1), path(fastq_2), val(project)
    val outdir

    output:
    tuple  val(project), path("*_MERGED.fastq.gz"), emit: merged_fastq

    """
    echo "Forward strand: $fastq_1"
    echo "Reverse strand: $fastq_2"
    echo "Project: $project"
    echo "Output path: $outdir/merged/$project"
    SeqPrep -f $fastq_1 -r $fastq_2 -1 1_unmerged.fastq.gz -2 2_unmerged.fastq.gz -s ${sampleId}_MERGED.fastq.gz
    """

}