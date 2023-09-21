
process EASEL {

    label 'light'
    container = '/hps/nobackup/rdf/metagenomics/singularity_cache_nextflow/quay.io-biocontainers-easel-0.48--hec16e2b_1.img'

    input:
    tuple val(meta), path(reads_fasta), path(cmsearch_deoverlap_out)

    output:
    tuple val(meta), path("*.easel_extracted.fasta"), emit: easel_coords
    tuple val(meta), path("*.matched_seqs_with_coords"), emit: matched_seqs_with_coords

    script:
    def fasta_unzipped = reads_fasta.name.replace(".gz", "")


    """
    gzip -c -d $reads_fasta > $fasta_unzipped
    awk '{print \$1"-"\$3"/"\$8"-"\$9" "\$8" "\$9" "\$1}' $cmsearch_deoverlap_out > ${meta.id}.matched_seqs_with_coords
    esl-sfetch --index $fasta_unzipped
    esl-sfetch -Cf $fasta_unzipped ${meta.id}.matched_seqs_with_coords > ${meta.id}.easel_extracted.fasta
    """

}
