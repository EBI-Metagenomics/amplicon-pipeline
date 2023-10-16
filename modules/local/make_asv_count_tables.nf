
process MAKE_ASV_COUNT_TABLES {

    label 'light'
    // publishDir "${outdir}/${project}/${sampleId}/asv-gen/${var_region}", mode : "copy" 

    input:
    tuple val(meta), val(var_region), path(maps), path(chimeric_prop), path(matched_prop), path(taxa), path(extracted_var_path), path(reads)

    output:
    tuple val(meta), val(var_region), path("*asv_krona_counts.txt"), emit: asv_count_tables_out

    """
    zcat $fastq_1 | sed -n "1~4p" > headers.txt
    python /hps/software/users/rdf/metagenomics/service-team/users/chrisata/asv_gen/bin/make_asv_count_table.py -t $taxa -f $map_1 -r $map_2 -a $extracted_var_path -hd ./headers.txt  -s $sampleId
    """

}
