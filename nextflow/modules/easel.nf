
process easel {

    label 'light'
    container = '/hps/nobackup/rdf/metagenomics/singularity_cache_nextflow/quay.io-biocontainers-easel-0.48--hec16e2b_1.img'

    input:
    tuple val(project), path(fasta), path(cmsearch_deoverlap_out)
    val outdir

    output:
    tuple val(project), path('*.matched_seqs_with_coords'), emit: easel_coords

    """
    awk '{print \$1"-"\$3"/"\$8"-"\$9" "\$8" "\$9" "\$1}' $cmsearch_deoverlap_out > ${fasta.baseName}.matched_seqs_with_coords
    """

}
