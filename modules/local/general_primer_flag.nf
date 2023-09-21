
process GENERAL_PRIMER_FLAG {
    // Check for the presence of primers in general
    
    label 'light'

    input:
    tuple val(meta), val(var_region), path(reads_merged)

    output:
    tuple val(meta), val(var_region), path("*general_primer_out.txt"), emit: general_primer_out

    """
    python /hps/software/users/rdf/metagenomics/service-team/users/chrisata/asv_gen/bin/are_there_primers_MERGED.py -i $reads_merged -s ${meta.id}_${var_region} -o ./
    """

}
