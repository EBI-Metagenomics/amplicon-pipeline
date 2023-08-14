nextflow.enable.dsl=2

// nextflow run nextflow/primer_trimming.nf --path /hps/nobackup/rdf/metagenomics/service-team/users/chrisata/asv_datasets/ERP123542/raw --project ERP123542 --outdir /hps/nobackup/rdf/metagenomics/service-team/users/chrisata/asv_nf_testing/merged
// nextflow run nextflow/primer_trimming.nf --path /hps/nobackup/rdf/metagenomics/service-team/users/chrisata/asv_datasets/ERP122862_subset --project ERP122862 --outdir /hps/nobackup/rdf/metagenomics/service-team/users/chrisata/asv_nf_testing/merged
// nextflow run nextflow/primer_trimming.nf -profile lsf --path /hps/nobackup/rdf/metagenomics/service-team/users/chrisata/asv_datasets/ERP123542/raw --project ERP123542 --outdir /hps/nobackup/rdf/metagenomics/service-team/users/chrisata/asv_nf_testing



include { CLASSIFY_VAR_REGIONS } from './modules/classify_var_regions.nf'
include { PARSE_VAR_CLASSIFICATION } from './modules/parse_var_classification'
include { CUTADAPT } from './modules/cutadapt.nf'
include { CONCAT_PRIMERS } from './modules/concat_primers.nf'

include { QC } from './subworkflows/qc_swf.nf'
include { CMSEARCH_SUBWF } from './subworkflows/cmsearch_swf.nf'
include { MAPSEQ_OTU_KRONA as MAPSEQ_OTU_KRONA_SSU} from './subworkflows/mapseq_otu_krona_swf.nf'
include { MAPSEQ_OTU_KRONA as MAPSEQ_OTU_KRONA_LSU} from './subworkflows/mapseq_otu_krona_swf.nf'
include { PRIMER_IDENTIFICATION } from './subworkflows/primer_identification_swf.nf'
include { AUTOMATIC_PRIMER_PREDICTION } from './subworkflows/automatic_primer_trimming.nf'


// Silva databases
// TODO: Move these to config
ssu_db_fasta = file("/hps/software/users/rdf/metagenomics/service-team/users/chrisata/asv_gen/data/silva_ssu-20200130/SSU.fasta")
ssu_db_tax = file("/hps/software/users/rdf/metagenomics/service-team/users/chrisata/asv_gen/data/silva_ssu-20200130/slv_ssu_filtered2.txt")
ssu_db_otu = file("/hps/software/users/rdf/metagenomics/service-team/users/chrisata/asv_gen/data/silva_ssu-20200130/ssu2.otu")
ssu_db_mscluster = file("/hps/software/users/rdf/metagenomics/service-team/users/chrisata/asv_gen/data/silva_ssu-20200130/SSU.fasta.mscluster")
ssu_label = "SSU"

lsu_db_fasta = file("/hps/software/users/rdf/metagenomics/service-team/users/chrisata/asv_gen/data/silva_lsu-20200130/LSU.fasta")
lsu_db_tax = file("/hps/software/users/rdf/metagenomics/service-team/users/chrisata/asv_gen/data/silva_lsu-20200130/slv_lsu_filtered2.txt")
lsu_db_otu = file("/hps/software/users/rdf/metagenomics/service-team/users/chrisata/asv_gen/data/silva_lsu-20200130/lsu2.otu")
lsu_db_mscluster = file("/hps/software/users/rdf/metagenomics/service-team/users/chrisata/asv_gen/data/silva_lsu-20200130/LSU.fasta.mscluster")
lsu_label = "LSU"

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

    PARSE_VAR_CLASSIFICATION(
        CLASSIFY_VAR_REGIONS.out.classify_var_summary,
        CLASSIFY_VAR_REGIONS.out.classify_var_regions,
        outdir
    )

    // Automatic primer identification
    PRIMER_IDENTIFICATION(
        QC.out.merged_reads,
        outdir
    )

    // Join flags with merged fastq files
    auto_trimming_input = PRIMER_IDENTIFICATION.out.conductor_out
                          .join(QC.out.merged_reads, by: [0, 1])


    // Run subworkflow for automatic primer trimming
    // Outputs either empty fasta file or fasta file containing predicted primers
    AUTOMATIC_PRIMER_PREDICTION(
        auto_trimming_input,
        outdir
    )

    // Join auto flags to std flags and generated a concatenated fasta file containing primers to trim off
    // This can contain any valid combination of stranded std/auto primers 
    concat_input = PRIMER_IDENTIFICATION.out.std_primer_out
                   .join(AUTOMATIC_PRIMER_PREDICTION.out.auto_primer_trimming_out, by: [0, 1])
   
    CONCAT_PRIMERS(
        concat_input,
        outdir
    )

    // Join concatenated primers to the fastp-cleaned paired reads files and run cutadapt on them
    cutadapt_input = CONCAT_PRIMERS.out.concat_primers_out
                     .join(QC.out.fastp_cleaned_fastq, by: [0, 1])
   
    CUTADAPT(
        cutadapt_input,
        outdir
    )

    // Just some logging for myself, will delete this eventually
    // final_out = fastp.out.cleaned_fastq
    // .join(CUTADAPT.out.cutadapt_out, remainder: true)
    // .map( { if (it[4] == null) { tuple(it[0], it[1], it[2], it[3]) } else { tuple(it[0], it[1], it[4], it[5]) }} )


}