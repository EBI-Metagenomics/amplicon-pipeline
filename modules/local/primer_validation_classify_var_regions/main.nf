
process PRIMER_VALIDATION_CLASSIFY_VAR_REGIONS {
    tag "$meta.id"
    label 'very_light'
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     "https://depot.galaxyproject.org/singularity/mgnify-pipelines-toolkit:${params.mpt_version}":
    //     "biocontainers/mgnify-pipelines-toolkit:${params.mpt_version}" }"
    container "oras://community.wave.seqera.io/library/pip_mgnify-pipelines-toolkit:64432eed2c673687"

    input:
    tuple val(meta), path(cmsearch_deoverlap_out), path(concat_primers_fasta)

    output:
    tuple val(meta), path("*primer_validation.tsv"), emit: primer_validation_out
    tuple val(meta), path("*primers.fasta")        , emit: validated_primers
    path "versions.yml"                            , emit: versions

    script:
    """
    // TODO remove absolute path obviously, this is temporary
    python /hps/software/users/rdf/metagenomics/service-team/users/chrisata/mgnify-pipelines-toolkit/mgnify_pipelines_toolkit/analysis/amplicon/primer_val_classification.py -i $cmsearch_deoverlap_out -f $concat_primers_fasta -s ${meta.id}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mgnify-pipelines-toolkit: ${params.mpt_version}
    END_VERSIONS
    """
}
