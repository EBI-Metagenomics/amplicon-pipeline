
process GENERAL_PRIMER_FLAG {
    // Check for the presence of primers in general
    
    label 'light'

    input:
    tuple  val(project), val(sampleId), val(var_region), path(fasta)
    val outdir

    output:
    tuple val(project), val(sampleId), val(var_region), path("*general_primer_out.txt"), emit: general_primer_out
    // path "*auto_primers.fasta", emit: auto_primer_seq

    """
    python /hps/software/users/rdf/metagenomics/service-team/users/chrisata/asv_gen/bin/are_there_primers_MERGED.py -i $fasta -s ${fasta.simpleName}_${var_region} -o ./
    """

}
