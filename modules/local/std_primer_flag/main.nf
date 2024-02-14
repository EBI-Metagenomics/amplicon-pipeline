
process STD_PRIMER_FLAG {
    // Check for presence of standard library of primers (stored in ./data/standard_primers)
    tag "$meta.id"
    label 'light'
    // publishDir "${outdir}/${project}/${sampleId}/primer-identification", mode : "copy" 

    input:
    tuple val(meta), path(reads_merged)

    output:
    tuple val(meta), path("*std_primers.fasta"), emit: std_primer_out
    path "*std_primer_out.txt"

    """
    python /hps/software/users/rdf/metagenomics/service-team/users/chrisata/asv_gen/bin/standard_primer_agrep.py -i $reads_merged -s ${meta.id}_${meta.var_region} -o ./
    """

}
