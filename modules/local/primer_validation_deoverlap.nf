
process PRIMER_VALIDATION_DEOVERLAP {

    label 'light' 
    container = '/hps/nobackup/rdf/metagenomics/singularity_cache_nextflow/quay.io-biocontainers-perl-5.22.2.1.img'

    input:
    tuple val(meta), path(cmsearch_out)
    path(rfam_clan)

    output:
    tuple val(meta), path("${cmsearch_out}.deoverlapped"), emit: cmsearch_deoverlap_out

    """
    perl /hps/software/users/rdf/metagenomics/service-team/users/chrisata/asv_gen/bin/cmsearch-deoverlap_primerval.pl -s --clanin $rfam_clan $cmsearch_out
    """

}
