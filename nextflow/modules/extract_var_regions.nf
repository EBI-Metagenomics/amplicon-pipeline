process EXTRACT_VAR_REGIONS {

    label 'light'
    container = '/hps/nobackup/rdf/metagenomics/singularity_cache/quay.io_biocontainers_seqtk:1.3.sif'

    input:
    tuple val(project), val(sampleId), path(var_region_path), path(fastq)
    val outdir

    output:
    tuple val(project), val(sampleId), val(var_region), path("*.fastq.gz"), emit: extracted_var_out
    
    script:
    var_region = "${var_region_path.baseName.split('\\.')[1,2].join('-')}"

    """
    seqtk subseq $fastq $var_region_path > ${fastq.baseName}_extracted.fastq
    gzip ${fastq.baseName}_extracted.fastq
    """

}

    // seqtk subseq $fastq $var_region_path > ${fastq.simpleName}_