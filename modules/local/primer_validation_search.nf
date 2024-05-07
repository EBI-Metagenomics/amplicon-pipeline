

process PRIMER_VALIDATION_SEARCH {

    label 'high_cpu_low_mem' // Will likely need to give this task more CPUs 
    container = '/hps/nobackup/rdf/metagenomics/service-team/singularity-cache/quay.io-biocontainers-infernal-1.1.4--pl5321hec16e2b_1.img'

    input:
    tuple val(meta), path(primer_fasta)
    path(rfam)

    output:
    tuple val(meta), path("*cmsearch_matches.tbl"), emit: cmsearch_out

    """
    cmsearch -g --noali -o /dev/null --tblout ${meta.id}.cmsearch_matches.tbl $rfam $primer_fasta
    """
//     cmsearch --hmmonly --hmmF1 1 --hmmF2 1 --hmmF3 1 -o /dev/null --tblout ${meta.id}.cmsearch_matches.tbl $rfam $primer_fasta

}
