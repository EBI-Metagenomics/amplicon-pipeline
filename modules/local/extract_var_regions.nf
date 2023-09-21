process EXTRACT_VAR_REGIONS {

    label 'light'
    // publishDir "${outdir}/${project}/${sampleId}/amplified-region-inference", pattern : "*.fastq.gz" , mode : "copy"
    container = '/hps/nobackup/rdf/metagenomics/singularity_cache/quay.io_biocontainers_seqtk:1.3.sif'

    input:
    tuple val(meta), path(var_region_path), path(reads_merged)

    output:
    tuple val(meta), val(var_region), path("*.fastq.gz"), emit: extracted_var_out
    tuple val(meta), val(var_region), path(var_region_path), emit: extracted_var_path
    
    script:
    var_region = "${var_region_path.baseName.split('\\.')[1,2].join('-')}"

    """
    seqtk subseq $reads_merged $var_region_path > ${meta.id}_${var_region}_extracted.fastq
    gzip ${meta.id}_${var_region}_extracted.fastq
    """

}