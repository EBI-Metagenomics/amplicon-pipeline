
process MAKE_ASV_COUNT_TABLES {

    label 'light'
    // publishDir "${outdir}/${project}", mode : "copy"

    input:
    tuple val(project), val(sampleId), val(var_region), path(map_1), path(map_2), path(chimeric_prop), path(matched_prop), path(taxa), path(extracted_var_path), path(fastq_1), path(fastq_2)
    val outdir

    output:
    tuple val(project), val(sampleId), val(var_region), path("*asv_krona_counts.txt"), emit: asv_count_tables_out

    """
    zcat $fastq_1 | sed -n "1~4p" > headers.txt
    python /hps/software/users/rdf/metagenomics/service-team/users/chrisata/asv_gen/bin/make_asv_count_table.py -t $taxa -f $map_1 -r $map_2 -a $extracted_var_path -hd ./headers.txt  -s $sampleId
    """

}
