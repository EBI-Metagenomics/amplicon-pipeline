
rfam = file("/hps/software/users/rdf/metagenomics/service-team/users/chrisata/asv_gen/data/rfam/ribo.cm")

process cmsearch {

    label 'light' // Will likely need to give this task more CPUs 
    container = '/hps/nobackup/rdf/metagenomics/singularity_cache_nextflow/quay.io-biocontainers-infernal-1.1.4--pl5321hec16e2b_1.img'

    input:
    tuple val(project), path(fasta)
    val outdir

    output:
    tuple val(project), path("*cmsearch_matches.tbl"), emit: cmsearch_out

    """
    cmsearch --cut_ga --noali --hmmonly -Z 1000 -o /dev/null --tblout ${fasta.simpleName}.cmsearch_matches.tbl $rfam $fasta
    """

}
