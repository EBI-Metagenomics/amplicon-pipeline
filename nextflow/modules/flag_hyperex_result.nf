// grep '^>' merged/ERP119447/ERR6093661_MERGED_hyperex.fa | cut -d' ' -f3-4 | sort | uniq -c
//  159923 forward=GTGCCAGCMGCCGCGGTAA reverse=GGACTACHVGGGTWTCTAAT
// grep -c '^>' merged/ERP119447/ERR6093661_MERGED_hyperex.fa 

// /hps/nobackup/rdf/metagenomics/service-team/users/chrisata/agrep/agrep -1 'GTGCCAGC[A,C]GCCGCGGTAA' merged/ERP111358/ERR2832485_MERGED.fasta

process flag_hyperex_result {

    // publishDir "${outdir}/merged/${project}", mode : "copy"

    input:
    tuple path(hyperex_fa), path(hyperex_gff), val(project)
    val outdir

    output: 
    tuple path("*_hyperex.fa"), path("*_hyperex.gff"), val(project), emit: hyperex_out

    """
    primer_count=$(grep '^>' $hyperex_fa | cut -d' ' -f3-4 | sort | uniq -c)
    read_count = $(grep -c '^>' merged/ERP119447/ERR6093661_MERGED_hyperex.fa)
    """
}
