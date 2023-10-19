
process DADA2 {
    // Run DADA2 pipeline including read-tracking

    label 'heavy'
    // publishDir "${outdir}/${project}/${sampleId}/asv-gen", pattern : "*.tsv" , mode : "copy" 
    // publishDir "${outdir}/${project}/${sampleId}/asv-gen", pattern : "*chimeric.txt" , mode : "copy" 
    // publishDir "${outdir}/${project}/${sampleId}/asv-gen", pattern : "*matched.txt" , mode : "copy" 


    input:
    tuple val(meta), val(var_region), path(reads)
    path silva_dada2_db

    output:
    tuple val(meta), val(var_region), path("*map.txt"), path("*chimeric.txt"), path("*matched.txt"), path("*taxa.tsv"), optional: true, emit: dada2_out

    """
    if [[ ${meta.single_end} = true ]]; then
        Rscript /hps/software/users/rdf/metagenomics/service-team/users/chrisata/asv_gen/bin/dada2.R ${meta.id} $silva_dada2_db $reads
    else
        Rscript /hps/software/users/rdf/metagenomics/service-team/users/chrisata/asv_gen/bin/dada2.R ${meta.id} $silva_dada2_db ${reads[0]} ${reads[1]}
    fi
    """

}
