
process GENERAL_PRIMER_FLAG {
    // Check for the presence of primers in general
    
    label 'light'

    input:
    tuple  val(project), val(sampleId), val(var_region), path(fastq)
    val outdir

    output:
    tuple val(project), val(sampleId), val(var_region), path("*general_primer_out.txt"), emit: general_primer_out

    """
    python /hps/software/users/rdf/metagenomics/service-team/users/chrisata/asv_gen/bin/are_there_primers_MERGED.py -i $fastq -s ${fastq.simpleName}_${var_region} -o ./
    """

}
