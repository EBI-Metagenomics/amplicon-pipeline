
process PRIMER_VALIDATION_CLASSIFY_VAR_REGIONS {

    label 'light' 

    input:
    tuple val(project), val(sampleId), path(cmsearch_deoverlap_out), path(concat_primers_fasta)
    val outdir

    output:
    tuple val(project), val(sampleId), path("*primer_validation.tsv"), emit: primer_validation_out

    """
    python /hps/software/users/rdf/metagenomics/service-team/users/chrisata/asv_gen/bin/primer_validation_classify_var_regions.py -i $cmsearch_deoverlap_out -f $concat_primers_fasta -s $sampleId 
    """

}
