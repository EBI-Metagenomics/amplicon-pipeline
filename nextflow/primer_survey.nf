nextflow.enable.dsl=2

// nextflow run nextflow/primer_survey.nf --path /hps/nobackup/rdf/metagenomics/service-team/users/chrisata/asv_nf_testing/data --outdir /hps/nobackup/rdf/metagenomics/service-team/users/chrisata/asv_nf_testing
// nextflow run nextflow/primer_survey.nf -profile lsf --path /hps/nobackup/rdf/metagenomics/service-team/users/chrisata/asv_datasets/survey_datasets/final_survey_dataset_data/data --outdir /hps/nobackup/rdf/metagenomics/service-team/users/chrisata/asv_datasets/survey_datasets/final_survey_dataset_data/merged

include { fastp } from './modules/fastp.nf'
include { seqprep_merge } from './modules/seqprep_merge.nf'
include { std_primer_flag } from './modules/std_primer_flag.nf'
include { general_primer_flag } from './modules/general_primer_flag.nf'


params.path = null
params.outdir = null
params.cpu = null

reads = Channel.fromFilePairs( "${params.path}/**/raw/*_{1,2}.fastq.gz" )
projects = Channel.fromPath( "${params.path}/*", type: 'dir' )
outdir = params.outdir
cpu = params.cpu

process COMBINE_READS_PROJECTS {

    label 'light'

    input:
    tuple val(sampleId), path(fastq)
    val proj_path

    output:
    tuple val(sampleId), path(fastq), val(project), emit: combined_reads_projects

    script:
    project = "${proj_path}".split('/').last()

    """
    """

}


workflow {
    

    COMBINE_READS_PROJECTS(reads, projects)
    fastp(COMBINE_READS_PROJECTS.out.combined_reads_projects, outdir)
    seqprep_merge(fastp.out.cleaned_fastq, outdir)
    std_primer_flag(seqprep_merge.out.merged_fastq, outdir)
    general_primer_flag(seqprep_merge.out.merged_fastq, outdir)


    
    
    // fastq_to_fasta(SeqPrep.out.merged_fastq, outdir)

    // std_primer_out = std_primer_out_parse(std_primer_agrep.out)
    
    // general_primer_out = general_primer_out_parse(general_primer_flag.out)

    // std_primer_out.view()
    // std_primer_out.count('2').view()
    // std_primer_out.count('1').view()
    // std_primer_out.count('0').view()
    
    // // general_primer_out.view()
    // general_primer_out.count('2').view()
    // general_primer_out.count('1').view()
    // general_primer_out.count('0').view()

}