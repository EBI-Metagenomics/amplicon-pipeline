
process general_primer_flag {

    label 'light'
    publishDir "${outdir}/merged/${project}", mode : "copy"
    // container = '/hps/nobackup/rdf/metagenomics/singularity_cache_nextflow/quay.io-biocontainers-biopython-1.75.img'

    input:
    tuple path(fasta), val(project)
    val outdir

    output:
    path "*general_primer_out.txt", emit: general_primer_out

    """
    python /hps/software/users/rdf/metagenomics/service-team/users/chrisata/asv_gen/bin/are_there_primers_MERGED.py -i $fasta -s ${fasta.simpleName} -o ./
    """

}
