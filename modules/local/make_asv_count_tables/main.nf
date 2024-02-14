
process MAKE_ASV_COUNT_TABLES {

    label 'light'
    // publishDir "${outdir}/${project}/${sampleId}/asv-gen/${var_region}", mode : "copy" 

    input:
    tuple val(meta), path(maps), path(taxa), path(extracted_var_path), path(reads)

    output:
    tuple val(meta), path("*asv_krona_counts.txt"), emit: asv_count_tables_out

    """
    if [[ ${meta.single_end} = true ]]; then
        zcat $reads | sed -n "1~4p" > headers.txt
        python /hps/software/users/rdf/metagenomics/service-team/users/chrisata/asv_gen/bin/make_asv_count_table.py -t $taxa -f $maps -a $extracted_var_path -hd ./headers.txt  -s ${meta.id}
    else
        zcat ${reads[0]} | sed -n "1~4p" > headers.txt
        python /hps/software/users/rdf/metagenomics/service-team/users/chrisata/asv_gen/bin/make_asv_count_table.py -t $taxa -f ${maps[0]} -r ${maps[1]} -a $extracted_var_path -hd ./headers.txt  -s ${meta.id}
    fi    
    """

}
