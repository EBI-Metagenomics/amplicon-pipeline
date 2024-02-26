
process PRIMER_VALIDATION_CLASSIFY_VAR_REGIONS {
    tag "$meta.id"
    label 'light' 

    input:
    tuple val(meta), path(cmsearch_deoverlap_out), path(concat_primers_fasta)

    output:
    tuple val(meta), path("*primer_validation.tsv"), emit: primer_validation_out

    """
    python /hps/software/users/rdf/metagenomics/service-team/users/chrisata/asv_gen/bin/primer_validation_classify_var_regions.py -i $cmsearch_deoverlap_out -f $concat_primers_fasta -s ${meta.id} 
    """

}
