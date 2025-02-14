
process PRIMER_VALIDATION_CLASSIFY_VAR_REGIONS {
    tag "$meta.id"
    label 'very_light'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "https://depot.galaxyproject.org/singularity/mgnify-pipelines-toolkit:${params.mpt_version}":
        "biocontainers/mgnify-pipelines-toolkit:${params.mpt_version}" }"

    input:
    tuple val(meta), path(cmsearch_deoverlap_out), path(concat_primers_fasta)

    output:
    tuple val(meta), path("*primer_validation.tsv"), emit: primer_validation_out
    tuple val(meta), path("*_primers.fasta")       , emit: final_concat_primers
    path "versions.yml"                            , emit: versions

    script:
    """
    primer_val_classification -i $cmsearch_deoverlap_out -f $concat_primers_fasta -s ${meta.id}

    num_primers=\$(grep -c '>' ${concat_primers_fasta})
    num_lines=\$(wc -l ${meta.id}_primer_validation.tsv | cut -d' ' -f1)
    num_lines=\$((\$num_lines - 1))

    cp ${concat_primers_fasta} ${meta.id}_primers.fasta

    if [ \$num_primers -lt 2 ] || [ \$num_lines -ne \$num_primers ]
    then
        echo "Primer validation didn't pass. Outputting empty file."
        rm ${meta.id}_primer_validation.tsv ${meta.id}_primers.fasta
        echo -n > ${meta.id}_primer_validation.tsv
        echo -n > ${meta.id}_primers.fasta
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mgnify-pipelines-toolkit: ${params.mpt_version}
    END_VERSIONS
    """
}
