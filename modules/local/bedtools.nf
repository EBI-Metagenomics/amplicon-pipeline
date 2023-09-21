
process BEDTOOLS {

    label 'light'
    container = '/hps/nobackup/rdf/metagenomics/singularity_cache_nextflow/quay.io-biocontainers-bedtools-2.23.0--h5b5514e_6.img'

    input:
    tuple val(meta), path(reads_fasta), path(maskfile)

    output:
    tuple val(meta), path("*ITS_masked.fasta"), emit: its_masked_out

    script:
    def fasta_unzipped = reads_fasta.name.replace(".gz", "")

    """
    gzip -c -d $reads_fasta > $fasta_unzipped
    bedtools maskfasta -fi $fasta_unzipped -bed $maskfile -fo ${meta.id}_ITS_masked.fasta
    """

}
