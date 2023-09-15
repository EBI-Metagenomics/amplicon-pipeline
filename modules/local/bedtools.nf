
process BEDTOOLS {

    label 'light'
    container = '/hps/nobackup/rdf/metagenomics/singularity_cache_nextflow/quay.io-biocontainers-bedtools-2.23.0--h5b5514e_6.img'


    input:
    tuple val(project), val(sampleId), path(fasta), path(maskfile)
    val outdir

    output:
    tuple val(project), val(sampleId), path("*ITS_masked.fasta"), emit: its_masked_out

    """
    echo "hellow"
    bedtools maskfasta -fi $fasta -bed $maskfile -fo ${sampleId}_ITS_masked.fasta
    """

}
