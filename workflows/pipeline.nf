
// nextflow run -resume main.nf --path /hps/nobackup/rdf/metagenomics/service-team/users/chrisata/asv_datasets/ERP122862_subset --project ERP122862 --outdir /hps/nobackup/rdf/metagenomics/service-team/users/chrisata/asv_nf_testing/merged
// nextflow run -resume main.nf --path samplesheet.csv --project ERP122862 --outdir /hps/nobackup/rdf/metagenomics/service

include { INPUT_CHECK } from '../subworkflows/local/input_check.nf'
include { READS_QC } from '../subworkflows/ebi-metagenomics/reads_qc/main.nf'
// include { CMSEARCH_SUBWF } from '../subworkflows/local/cmsearch_swf.nf'
// include { ITS_SWF } from '../subworkflows/local/its_swf.nf'
// include { MAPSEQ_OTU_KRONA as MAPSEQ_OTU_KRONA_SSU} from '../subworkflows/local/mapseq_otu_krona_swf.nf'
// include { MAPSEQ_OTU_KRONA as MAPSEQ_OTU_KRONA_LSU} from '../subworkflows/local/mapseq_otu_krona_swf.nf'
// include { MAPSEQ_OTU_KRONA as MAPSEQ_OTU_KRONA_UNITE} from '../subworkflows/local/mapseq_otu_krona_swf.nf'
// include { MAPSEQ_OTU_KRONA as MAPSEQ_OTU_KRONA_ITSONEDB} from '../subworkflows/local/mapseq_otu_krona_swf.nf'
// include { AMP_REGION_INFERENCE } from '../subworkflows/local/amp_region_inference_swf.nf'
// include { PRIMER_IDENTIFICATION } from '../subworkflows/local/primer_identification_swf.nf'
// include { PRIMER_VALIDATION } from '../subworkflows/local/primer_validation_swf.nf'
// include { AUTOMATIC_PRIMER_PREDICTION } from '../subworkflows/local/automatic_primer_prediction.nf'
// include { CONCAT_PRIMER_CUTADAPT } from '../subworkflows/local/concat_primer_cutadapt.nf'
// include { DADA2_KRONA } from '../subworkflows/local/dada2_krona_swf.nf'

// Initialise different database inputs for MapSeq+Krona
ssu_mapseq_krona_tuple = tuple(file(params.ssu_db_fasta), file(params.ssu_db_tax), file(params.ssu_db_otu), file(params.ssu_db_mscluster), params.ssu_label)
lsu_mapseq_krona_tuple = tuple(file(params.lsu_db_fasta), file(params.lsu_db_tax), file(params.lsu_db_otu), file(params.lsu_db_mscluster), params.lsu_label)
itsonedb_mapseq_krona_tuple = tuple(file(params.itsone_db_fasta), file(params.itsone_db_tax), file(params.itsone_db_otu), file(params.itsone_db_mscluster), params.itsone_label)
unite_mapseq_krona_tuple = tuple(file(params.unite_db_fasta), file(params.unite_db_tax), file(params.unite_db_otu), file(params.unite_db_mscluster), params.unite_label)

// Initialise database inputs for DADA2+Krona
silva_dada2_db = file(params.silva_dada2_db)
dada2_krona_tuple = tuple(file(params.ssu_db_fasta), file(params.ssu_db_tax), file(params.ssu_db_otu), file(params.ssu_db_mscluster), params.dada2_silva_label)

params.path = null
params.project = null
params.outdir = null
params.mode = null


project = Channel.value( params.project )


samplesheet = file( params.path )
// reads = Channel.fromFilePairs( "${params.path}/*{1,2}.fastq.gz", size: 2)
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


workflow AMPLICON_PIPELINE_V6 {

    // TODO: need to get this working for single-end reads
    // TODO: need to change the 'project' id in the tuples to meta to work with nf-core
    // TODO: investigate primer val deoverlap script more


    INPUT_CHECK(samplesheet)

    // Quality control

    READS_QC(
        INPUT_CHECK.out.reads
    )
    
    READS_QC.out.reads_se_and_merged.view()

    // // Cmsearch subworkflow to find rRNA reads for SSU+LSU
    // CMSEARCH_SUBWF(
    //     QC.out.merged_fasta,
    //     outdir
    // )

    // // Masking subworkflow to find rRNA reads for ITS
    // ITS_SWF(
    //     QC.out.merged_fasta,
    //     CMSEARCH_SUBWF.out.concat_ssu_lsu_coords,
    //     outdir
    // )

    // // Next four subworkflow calls are MapSeq annotation + Krona generation for SSU+LSU+ITS
    // MAPSEQ_OTU_KRONA_SSU(
    //     CMSEARCH_SUBWF.out.ssu_fasta,
    //     ssu_mapseq_krona_tuple,
    //     outdir
    // )

    // MAPSEQ_OTU_KRONA_LSU(
    //     CMSEARCH_SUBWF.out.lsu_fasta,
    //     lsu_mapseq_krona_tuple,
    //     outdir
    // )    

    // MAPSEQ_OTU_KRONA_ITSONEDB(
    //     ITS_SWF.out.its_masked_out,
    //     itsonedb_mapseq_krona_tuple,
    //     outdir
    // )    

    // MAPSEQ_OTU_KRONA_UNITE(
    //     ITS_SWF.out.its_masked_out,
    //     unite_mapseq_krona_tuple,
    //     outdir
    // )    

    // // Infer amplified variable regions for SSU, extract reads for each amplified region if there are more than one
    // AMP_REGION_INFERENCE(
    //     CMSEARCH_SUBWF.out.cmsearch_deoverlap_out,
    //     QC.out.merged_reads,
    //     outdir
    // )

    // // Identify whether primers exist or not in reads, separated by different amplified regions if more than one exists in a run
    // PRIMER_IDENTIFICATION(
    //     AMP_REGION_INFERENCE.out.extracted_var_out,
    //     outdir
    // )

    // // Join primer identification flags with reads belonging to each run+amp_region
    // auto_trimming_input = PRIMER_IDENTIFICATION.out.conductor_out
    //                       .join(AMP_REGION_INFERENCE.out.extracted_var_out, by: [0, 1, 2])


    // // Run subworkflow for automatic primer prediction
    // // Outputs empty fasta file if no primers, or fasta file containing predicted primers
    // AUTOMATIC_PRIMER_PREDICTION(
    //     auto_trimming_input,
    //     outdir
    // )

    // // Concatenate the different combinations of stranded std/auto primers for each run+amp_region
    // concat_input = PRIMER_IDENTIFICATION.out.std_primer_out
    //                .join(AUTOMATIC_PRIMER_PREDICTION.out.auto_primer_trimming_out, by: [0, 1, 2])

    // // Concatenate all primers for for a run, send them to cutadapt with original QCd reads for primer trimming
    // CONCAT_PRIMER_CUTADAPT(
    //     concat_input,
    //     QC.out.fastp_cleaned_fastq,
    //     outdir
    // )

    // primer_validation_input = CONCAT_PRIMER_CUTADAPT.out.final_concat_primers_out
    //                           .map{ tuple(it[0], it[1], it[3]) }

    // // Verify that any identified primers (both std+auto) actually match to regions of the SSU gene (for Bacteria/Archaea/Eukaryotes)
    // // Output of this (a .tsv file) will go to CDCH
    // PRIMER_VALIDATION(
    //     primer_validation_input,
    //     outdir
    // )

    // // Prepare DADA2 input (either fastp reads if no primer trimming was done, or cutadapt output if primers were trimmed)
    // dada2_input = CONCAT_PRIMER_CUTADAPT.out.cutadapt_out
    //               .join(QC.out.fastp_cleaned_fastq, by: [0, 1])
    //               .map( { if (it[2] != null && it[3] == null) { tuple(it[0], it[1], it[2], it[5], it[6]) } else { tuple(it[0], it[1], it[2], it[3], it[4]) }} )

    // // Run DADA2 ASV generation + generate Krona plots for each run+amp_region 
    // DADA2_KRONA(
    //     dada2_input,
    //     AMP_REGION_INFERENCE.out.concat_var_regions,
    //     AMP_REGION_INFERENCE.out.extracted_var_path,
    //     QC.out.fastp_cleaned_fastq,
    //     silva_dada2_db,
    //     dada2_krona_tuple,
    //     outdir
    // )

}