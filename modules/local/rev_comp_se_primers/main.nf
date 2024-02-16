
process REV_COMP_SE_PRIMERS {
    tag "$meta.id"
    label 'light'
    // publishDir "${outdir}/${project}/${sampleId}/primer-identification", mode : "copy" 
    
    input:
    tuple val(meta), path(final_concat_primers)

    output:
    tuple val(meta), path("*rev_comp_se_primers.fasta"), optional:true, emit: rev_comp_se_primers_out

    """
    if [[ -s $final_concat_primers && ${meta.single_end} = true ]]; then
        python /hps/software/users/rdf/metagenomics/service-team/users/chrisata/asv_gen/bin/rev_comp_se_primers.py -i $final_concat_primers -s ${meta.id} -o ./
    elif [[ -s $final_concat_primers && ${meta.single_end} = false ]]; then
        cat $final_concat_primers > ./${meta.id}_rev_comp_se_primers.fasta
    else
        touch ./${meta.id}_rev_comp_se_primers.fasta
    fi
    """

}
