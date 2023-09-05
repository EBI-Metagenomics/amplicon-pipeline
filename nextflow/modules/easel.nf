
process EASEL {

    label 'light'
    container = '/hps/nobackup/rdf/metagenomics/singularity_cache_nextflow/quay.io-biocontainers-easel-0.48--hec16e2b_1.img'

    input:
    tuple val(project), val(sampleId), path(fasta), path(cmsearch_deoverlap_out)
    val outdir

    output:
    tuple val(project), val(sampleId), path("*.easel_extracted.fasta"), emit: easel_coords
    path("*.matched_seqs_with_coords"), emit: matched_seqs_with_coords

    """
    awk '{print \$1"-"\$3"/"\$8"-"\$9" "\$8" "\$9" "\$1}' $cmsearch_deoverlap_out > ${fasta.baseName}.matched_seqs_with_coords
    esl-sfetch --index $fasta
    esl-sfetch -Cf $fasta ${fasta.baseName}.matched_seqs_with_coords > ${fasta.baseName}.easel_extracted.fasta
    """

}
