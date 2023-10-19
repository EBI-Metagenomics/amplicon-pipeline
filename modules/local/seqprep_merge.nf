
process SEQPREP_MERGE {
    // Merge paired-end reads with seqprep

    label 'light'
    container = '/hps/nobackup/rdf/metagenomics/service-team/singularity-cache/quay.io-biocontainers-seqprep-1.3.2--hed695b0_4.img'

    input:
    tuple val(project), val(sampleId), path(fastq_1), path(fastq_2)
    val outdir

    output:
    tuple val(project), val(sampleId), path("*_MERGED.fastq.gz"), emit: merged_fastq

    """
    SeqPrep -f $fastq_1 -r $fastq_2 -1 1_unmerged.fastq.gz -2 2_unmerged.fastq.gz -s ${sampleId}_MERGED.fastq.gz
    """

}