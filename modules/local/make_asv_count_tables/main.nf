
process MAKE_ASV_COUNT_TABLES {
    tag "$meta.id"
    label 'light'
    conda "bioconda::mgnify-pipelines-toolkit=0.1.0"

    input:
    tuple val(meta), path(maps), path(taxa), path(extracted_var_path), path(reads)

    output:
    tuple val(meta), path("*asv_krona_counts.txt"), emit: asv_count_tables_out

    """
    if [[ ${meta.single_end} = true ]]; then
        zcat $reads | sed -n "1~4p" > headers.txt
        make_asv_count_table -t $taxa -f $maps -a $extracted_var_path -hd ./headers.txt  -s ${meta.id}
    else
        zcat ${reads[0]} | sed -n "1~4p" > headers.txt
        make_asv_count_table -t $taxa -f ${maps[0]} -r ${maps[1]} -a $extracted_var_path -hd ./headers.txt  -s ${meta.id}
    fi    
    """

}
