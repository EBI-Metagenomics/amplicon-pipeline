nextflow.enable.dsl=2

// nextflow run nextflow/primer_trimming.nf --path /hps/nobackup/rdf/metagenomics/service-team/users/chrisata/asv_datasets/ERP123542/raw --project ERP123542 --outdir /hps/nobackup/rdf/metagenomics/service-team/users/chrisata/asv_nf_testing/merged
// nextflow run nextflow/primer_trimming.nf --path /hps/nobackup/rdf/metagenomics/service-team/users/chrisata/asv_datasets/ERP122862_subset --project ERP122862 --outdir /hps/nobackup/rdf/metagenomics/service-team/users/chrisata/asv_nf_testing/merged
// nextflow run nextflow/primer_trimming.nf -profile lsf --path /hps/nobackup/rdf/metagenomics/service-team/users/chrisata/asv_datasets/ERP123542/raw --project ERP123542 --outdir /hps/nobackup/rdf/metagenomics/service-team/users/chrisata/asv_nf_testing



include { CLASSIFY_VAR_REGIONS } from './modules/classify_var_regions.nf'
include { PARSE_VAR_CLASSIFICATION } from './modules/parse_var_classification.nf'
include { EXTRACT_VAR_REGIONS } from './modules/extract_var_regions.nf'
include { CUTADAPT } from './modules/cutadapt.nf'
include { CONCAT_PRIMERS } from './modules/concat_primers.nf'
include { FINAL_CONCAT_PRIMERS } from './modules/final_concat_primers.nf'
include { DADA2 } from './modules/dada2.nf'
include { MAKE_ASV_COUNT_TABLES } from './modules/make_asv_count_tables.nf'
include { KRONA } from './modules/krona.nf'


include { QC } from './subworkflows/qc_swf.nf'
include { CMSEARCH_SUBWF } from './subworkflows/cmsearch_swf.nf'
include { MAPSEQ_OTU_KRONA as MAPSEQ_OTU_KRONA_SSU} from './subworkflows/mapseq_otu_krona_swf.nf'
include { MAPSEQ_OTU_KRONA as MAPSEQ_OTU_KRONA_LSU} from './subworkflows/mapseq_otu_krona_swf.nf'
include { PRIMER_IDENTIFICATION } from './subworkflows/primer_identification_swf.nf'
include { AUTOMATIC_PRIMER_PREDICTION } from './subworkflows/automatic_primer_trimming.nf'


// TODO: Move these to config
// Silva databases
ssu_db_fasta = file("/hps/nobackup/rdf/metagenomics/service-team/users/chrisata/silva_ssu-20200130/SSU.fasta")
ssu_db_tax = file("/hps/nobackup/rdf/metagenomics/service-team/users/chrisata/silva_ssu-20200130/slv_ssu_filtered2.txt")
ssu_db_otu = file("/hps/nobackup/rdf/metagenomics/service-team/users/chrisata/silva_ssu-20200130/ssu2.otu")
ssu_db_mscluster = file("/hps/nobackup/rdf/metagenomics/service-team/users/chrisata/silva_ssu-20200130/SSU.fasta.mscluster")
ssu_label = "SSU"

lsu_db_fasta = file("/hps/nobackup/rdf/metagenomics/service-team/users/chrisata/silva_lsu-20200130/LSU.fasta")
lsu_db_tax = file("/hps/nobackup/rdf/metagenomics/service-team/users/chrisata/silva_lsu-20200130/slv_lsu_filtered2.txt")
lsu_db_otu = file("/hps/nobackup/rdf/metagenomics/service-team/users/chrisata/silva_lsu-20200130/lsu2.otu")
lsu_db_mscluster = file("/hps/nobackup/rdf/metagenomics/service-team/users/chrisata/silva_lsu-20200130/LSU.fasta.mscluster")
lsu_label = "LSU"
silva_dada2_db = file("/hps/software/users/rdf/metagenomics/service-team/users/chrisata/asv_gen/data/silva_nr99_v138.1_train_set.fa.gz")

params.path = null
params.project = null
params.outdir = null

reads = Channel.fromFilePairs( "${params.path}/*_{1,2}.fastq.gz" )
project = Channel.value( params.project )
outdir = params.outdir

workflow {

    // Workflow for running the primer-trimming survey on sequencing reads

    // Quality control
    // TODO: need to get this working for single-end reads
    // TODO: need to change the 'project' id in the tuples to meta to work with nf-core
    // TODO: organise the publishDirs in a sensible way
    // TODO: rewrite arparse descriptions

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

    ssu_mapseq_krona_tuple = tuple(ssu_db_fasta, ssu_db_tax, ssu_db_otu, ssu_db_mscluster, ssu_label)
    lsu_mapseq_krona_tuple = tuple(lsu_db_fasta, lsu_db_tax, lsu_db_otu, lsu_db_mscluster, lsu_label)

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

    // Classify amplified regions
    CLASSIFY_VAR_REGIONS(
        CMSEARCH_SUBWF.out.cmsearch_deoverlap_out,
        outdir
    )

    extract_var_input = CLASSIFY_VAR_REGIONS.out.classify_var_regions
    .transpose()
    .combine(QC.out.merged_reads, by: [0, 1])

    EXTRACT_VAR_REGIONS(
        extract_var_input,
        outdir
    )


    // Automatic primer identification
    PRIMER_IDENTIFICATION(
        EXTRACT_VAR_REGIONS.out.extracted_var_out,
        outdir
    )

    // Join flags with merged fastq files
    auto_trimming_input = PRIMER_IDENTIFICATION.out.conductor_out
                          .join(EXTRACT_VAR_REGIONS.out.extracted_var_out, by: [0, 1, 2])


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
   
    CONCAT_PRIMERS(
        concat_input,
        outdir
    )

    final_concat_primers_input = CONCAT_PRIMERS.out.concat_primers_out
    .groupTuple(by: [0, 1])


    FINAL_CONCAT_PRIMERS(
        final_concat_primers_input,
        outdir
    )

    // Join concatenated primers to the fastp-cleaned paired reads files and run cutadapt on them
    cutadapt_input = FINAL_CONCAT_PRIMERS.out.final_concat_primers_out
                     .join(QC.out.fastp_cleaned_fastq, by: [0, 1])

    CUTADAPT(
        cutadapt_input,
        outdir
    )

    // Prepare DADA2 input (either fastp reads or cutadapt reads)

    // TODO: maybe put back the remainder=true for this join?
    dada2_input = CUTADAPT.out.cutadapt_out
                  .join(QC.out.fastp_cleaned_fastq, by: [0, 1])
                  .map( { if (it[2] != null && it[3] == null) { tuple(it[0], it[1], it[2], it[5], it[6]) } else { tuple(it[0], it[1], it[2], it[3], it[4]) }} )

    DADA2(
        dada2_input,
        silva_dada2_db,
        outdir
    )

    split_input = DADA2.out.dada2_out
                  .transpose()
                  .join(EXTRACT_VAR_REGIONS.out.extracted_var_path, by: [0, 1, 2])
                  

    multi_region_concats = split_input
    .join(CLASSIFY_VAR_REGIONS.out.concat_var_regions, by: [0, 1])
    .map( {tuple(it[0], it[1], "concat", it[3], it[4], it[5], it[6], it[7], it[10])} )
    
    
    final_asv_count_table_input = split_input
                                  .mix(multi_region_concats)
                                  .combine(QC.out.fastp_cleaned_fastq, by: [0, 1])

                                
    MAKE_ASV_COUNT_TABLES(
        final_asv_count_table_input,
        outdir
    )

    asv_krona_input = MAKE_ASV_COUNT_TABLES.out.asv_count_tables_out
                      .map( {it[0, 1, 3]} )

    KRONA(
        asv_krona_input,
        ssu_mapseq_krona_tuple,
        outdir
    )
    

}