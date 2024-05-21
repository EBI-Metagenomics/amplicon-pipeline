
process PRIMER_VALIDATION_CLASSIFY_VAR_REGIONS {
    tag "$meta.id"
    label 'very_light'
    conda "${projectDir}/conf/environment.yml"
    // TODO: uncomment container when you release fix to mpt
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     "https://depot.galaxyproject.org/singularity/mgnify-pipelines-toolkit:${params.mpt_version}":
    //     "biocontainers/mgnify-pipelines-toolkit:${params.mpt_version} }"

    input:
    tuple val(meta), path(cmsearch_deoverlap_out), path(concat_primers_fasta)

    output:
    tuple val(meta), path("*primer_validation.tsv"), emit: primer_validation_out

    script:
    """
    primer_validation_classify_var_regions.py -i $cmsearch_deoverlap_out -f $concat_primers_fasta -s ${meta.id} 
    """
}
