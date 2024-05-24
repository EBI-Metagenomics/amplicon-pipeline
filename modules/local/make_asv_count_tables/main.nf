
process MAKE_ASV_COUNT_TABLES {
    tag "$meta.id"
    label 'process_long'
    // TODO: uncomment container when you release fix to mpt
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     "https://depot.galaxyproject.org/singularity/mgnify-pipelines-toolkit:${params.mpt_version}":
    //     "biocontainers/mgnify-pipelines-toolkit:${params.mpt_version}" }"

    input:
    tuple val(meta), path(maps), path(asvtaxtable), path(reads), path(extracted_var_path)

    output:
    tuple val(meta), path("*asv_krona_counts.txt"), emit: asv_count_tables_out
    path "versions.yml"                           , emit: versions

    script:
    """
    if [[ ${meta.single_end} = true ]]; then
        zcat $reads | sed -n "1~4p" > headers.txt
        make_asv_count_table.py -t $asvtaxtable -f $maps -a $extracted_var_path -hd ./headers.txt  -s ${meta.id}
    else
        zcat ${reads[0]} | sed -n "1~4p" > headers.txt
        make_asv_count_table.py -t $asvtaxtable -f ${maps[0]} -r ${maps[1]} -a $extracted_var_path -hd ./headers.txt  -s ${meta.id}
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mgnify-pipelines-toolkit: ${params.mpt_version}
    END_VERSIONS
    """
}
