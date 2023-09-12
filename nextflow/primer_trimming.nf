nextflow.enable.dsl=2

// nextflow run nextflow/primer_trimming.nf --path /hps/nobackup/rdf/metagenomics/service-team/users/chrisata/asv_datasets/ERP123542/raw --project ERP123542 --outdir /hps/nobackup/rdf/metagenomics/service-team/users/chrisata/asv_nf_testing/merged
// nextflow run nextflow/primer_trimming.nf --path /hps/nobackup/rdf/metagenomics/service-team/users/chrisata/asv_datasets/ERP122862_subset --project ERP122862 --outdir /hps/nobackup/rdf/metagenomics/service-team/users/chrisata/asv_nf_testing/merged
// nextflow run nextflow/primer_trimming.nf -profile lsf --path /hps/nobackup/rdf/metagenomics/service-team/users/chrisata/asv_datasets/ERP123542/raw --project ERP123542 --outdir /hps/nobackup/rdf/metagenomics/service-team/users/chrisata/asv_nf_testing

include { QC } from './subworkflows/qc_swf.nf'
include { CMSEARCH_SUBWF } from './subworkflows/cmsearch_swf.nf'
include { ITS_SWF } from './subworkflows/its_swf.nf'
include { MAPSEQ_OTU_KRONA as MAPSEQ_OTU_KRONA_SSU} from './subworkflows/mapseq_otu_krona_swf.nf'
include { MAPSEQ_OTU_KRONA as MAPSEQ_OTU_KRONA_LSU} from './subworkflows/mapseq_otu_krona_swf.nf'
include { MAPSEQ_OTU_KRONA as MAPSEQ_OTU_KRONA_UNITE} from './subworkflows/mapseq_otu_krona_swf.nf'
include { MAPSEQ_OTU_KRONA as MAPSEQ_OTU_KRONA_ITSONEDB} from './subworkflows/mapseq_otu_krona_swf.nf'
include { AMP_REGION_INFERENCE } from './subworkflows/amp_region_inference_swf.nf'
include { PRIMER_IDENTIFICATION } from './subworkflows/primer_identification_swf.nf'
include { PRIMER_VALIDATION } from './subworkflows/primer_validation_swf.nf'
include { AUTOMATIC_PRIMER_PREDICTION } from './subworkflows/automatic_primer_trimming.nf'
include { CONCAT_PRIMER_CUTADAPT } from './subworkflows/concat_primer_cutadapt.nf'
include { DADA2_KRONA } from './subworkflows/dada2_krona_swf.nf'

ssu_mapseq_krona_tuple = tuple(file(params.ssu_db_fasta), file(params.ssu_db_tax), file(params.ssu_db_otu), file(params.ssu_db_mscluster), params.ssu_label)
lsu_mapseq_krona_tuple = tuple(file(params.lsu_db_fasta), file(params.lsu_db_tax), file(params.lsu_db_otu), file(params.lsu_db_mscluster), params.lsu_label)
itsonedb_mapseq_krona_tuple = tuple(file(params.itsone_db_fasta), file(params.itsone_db_tax), file(params.itsone_db_otu), file(params.itsone_db_mscluster), params.itsone_label)
unite_mapseq_krona_tuple = tuple(file(params.unite_db_fasta), file(params.unite_db_tax), file(params.unite_db_otu), file(params.unite_db_mscluster), params.unite_label)

silva_dada2_db = file(params.silva_dada2_db)

params.path = null
params.project = null
params.outdir = null
params.mode = null


project = Channel.value( params.project )
reads = Channel.fromFilePairs( "${params.path}/*{1,2}.fastq.gz", size: 2)
// single_reads = Channel.fromPath("${params.path}/*.fastq.gz")
outdir = params.outdir

// paired_read_paths = reads
// .map { it[1] }
// .collect()

// paired_read_paths.view()

// single_reads
// .map{ it in paired_read_paths }
// .view()
// .join(paired_read_paths)
// .view()


workflow {

    // Workflow for running the primer-trimming survey on sequencing reads

    // Quality control
    // TODO: need to get this working for single-end reads
    // TODO: need to change the 'project' id in the tuples to meta to work with nf-core
    // TODO: organise the publishDirs in a sensible way
    // TODO: rewrite arparse descriptions
    // TODO: investigate primer val deoverlap script more
    

    QC(
        project,
        reads,
        outdir
    )

    // cmsearch to get rRNA matches
    CMSEARCH_SUBWF(
        QC.out.merged_fasta,
        outdir
    )

    ITS_SWF(
        QC.out.merged_fasta,
        CMSEARCH_SUBWF.out.concat_ssu_lsu_coords,
        outdir
    )

    MAPSEQ_OTU_KRONA_SSU(
        CMSEARCH_SUBWF.out.ssu_fasta,
        ssu_mapseq_krona_tuple,
        outdir
    )

    MAPSEQ_OTU_KRONA_LSU(
        CMSEARCH_SUBWF.out.lsu_fasta,
        lsu_mapseq_krona_tuple,
        outdir
    )    

    MAPSEQ_OTU_KRONA_ITSONEDB(
        ITS_SWF.out.its_masked_out,
        itsonedb_mapseq_krona_tuple,
        outdir
    )    

    MAPSEQ_OTU_KRONA_UNITE(
        ITS_SWF.out.its_masked_out,
        unite_mapseq_krona_tuple,
        outdir
    )    

    // Classify amplified regions

    AMP_REGION_INFERENCE(
        CMSEARCH_SUBWF.out.cmsearch_deoverlap_out,
        QC.out.merged_reads,
        outdir
    )

    // Automatic primer identification
    PRIMER_IDENTIFICATION(
        AMP_REGION_INFERENCE.out.extracted_var_out,
        outdir
    )

    // Join flags with merged fastq files
    auto_trimming_input = PRIMER_IDENTIFICATION.out.conductor_out
                          .join(AMP_REGION_INFERENCE.out.extracted_var_out, by: [0, 1, 2])


    // Run subworkflow for automatic primer trimming
    // Outputs either empty fasta file or fasta file containing predicted primers
    AUTOMATIC_PRIMER_PREDICTION(
        auto_trimming_input,
        outdir
    )

    // Join auto flags to std flags and generated a concatenated fasta file containing primers to trim off
    // This can contain any valid combination of stranded std/auto primers 
    concat_input = PRIMER_IDENTIFICATION.out.std_primer_out
                   .join(AUTOMATIC_PRIMER_PREDICTION.out.auto_primer_trimming_out, by: [0, 1, 2])

    CONCAT_PRIMER_CUTADAPT(
        concat_input,
        QC.out.fastp_cleaned_fastq,
        outdir
    )

    // Prepare DADA2 input (either fastp reads or cutadapt reads)
    // TODO: maybe put back the remainder=true for this join?
    dada2_input = CONCAT_PRIMER_CUTADAPT.out.cutadapt_out
                  .join(QC.out.fastp_cleaned_fastq, by: [0, 1])
                  .map( { if (it[2] != null && it[3] == null) { tuple(it[0], it[1], it[2], it[5], it[6]) } else { tuple(it[0], it[1], it[2], it[3], it[4]) }} )

    DADA2_KRONA(
        dada2_input,
        AMP_REGION_INFERENCE.out.concat_var_regions,
        AMP_REGION_INFERENCE.out.extracted_var_path,
        QC.out.fastp_cleaned_fastq,
        silva_dada2_db,
        ssu_mapseq_krona_tuple,
        outdir
    )

    primer_validation_input = CONCAT_PRIMER_CUTADAPT.out.final_concat_primers_out
                              .map{ tuple(it[0], it[1], it[3]) }

    PRIMER_VALIDATION(
        primer_validation_input,
        outdir
    )

}