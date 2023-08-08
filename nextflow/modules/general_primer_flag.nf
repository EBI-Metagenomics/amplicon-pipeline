
process general_primer_flag {
    // Check for the presence of primers in general
    
    label 'light'
    publishDir "${outdir}/${project}", mode : "copy"
    // container = '/hps/nobackup/rdf/metagenomics/singularity_cache_nextflow/quay.io-biocontainers-biopython-1.75.img'

    input:
    tuple  val(project), path(fasta)
    val outdir

    output:
    tuple val(project), path("*general_primer_out.txt"), emit: general_primer_out
    // path "*auto_primers.fasta", emit: auto_primer_seq

    """
    python /hps/software/users/rdf/metagenomics/service-team/users/chrisata/asv_gen/bin/are_there_primers_MERGED.py -i $fasta -s ${fasta.simpleName} -o ./
    """

}
