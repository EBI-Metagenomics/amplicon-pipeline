
process DADA2 {
    // Run DADA2 pipeline including read-tracking

    label 'light'
    publishDir "${outdir}/${project}/${sampleId}/asv-gen", pattern : "*.tsv" , mode : "copy" 
    publishDir "${outdir}/${project}/${sampleId}/asv-gen", pattern : "*chimeric.txt" , mode : "copy" 
    publishDir "${outdir}/${project}/${sampleId}/asv-gen", pattern : "*matched.txt" , mode : "copy" 


    input:
    tuple val(project), val(sampleId), val(var_region), path(fastq_1), path(fastq_2)
    path silva_dada2_db
    val outdir

    output:
    tuple val(project), val(sampleId), val(var_region), path("*1_map.txt"), path("*2_map.txt"), path("*chimeric.txt"), path("*matched.txt"), path("*taxa.tsv"), emit: dada2_out

    """
    Rscript /hps/software/users/rdf/metagenomics/service-team/users/chrisata/asv_gen/bin/dada2.R $fastq_1 $fastq_2 $sampleId $silva_dada2_db
    """

}
