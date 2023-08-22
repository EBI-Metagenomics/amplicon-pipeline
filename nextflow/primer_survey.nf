nextflow.enable.dsl=2

// nextflow run nextflow/primer_survey.nf --path /hps/nobackup/rdf/metagenomics/service-team/users/chrisata/asv_nf_testing/data --outdir /hps/nobackup/rdf/metagenomics/service-team/users/chrisata/asv_nf_testing
// nextflow run nextflow/primer_survey.nf -profile lsf --path /hps/nobackup/rdf/metagenomics/service-team/users/chrisata/asv_datasets/survey_datasets/final_survey_dataset_data/data --outdir /hps/nobackup/rdf/metagenomics/service-team/users/chrisata/asv_datasets/survey_datasets/final_survey_dataset_data/merged

// include { fastp } from './modules/fastp.nf'
// include { seqprep_merge } from './modules/seqprep_merge.nf'
include { STD_PRIMER_FLAG } from './modules/std_primer_flag.nf'
// include { general_primer_flag } from './modules/general_primer_flag.nf'
// include { trimming_conductor } from './modules/trimming_conductor.nf'
// include { parse_conductor } from './modules/parse_conductor.nf'
// include { cutadapt } from './modules/cutadapt.nf'
// include { concat_primers } from './modules/concat_primers.nf'
// include { AUTOMATIC_PRIMER_TRIMMING } from './subworkflows/automatic_primer_trimming.nf'

params.path = null
params.outdir = null
params.filter = file("$projectDir/assets/NO_FILE")

reads = Channel.fromPath( "${params.path}/**/*MERGED.fastq.gz" )
projects = Channel.fromPath( "${params.path}/*", type: 'dir' )
outdir = params.outdir

process COMBINE_READS_PROJECTS {
    // Process for combining Project ID and paired-end reads into one tuple

    label 'light'

    input:
    // tuple val(sampleId), path(fastq)
    path fastq
    val proj_path

    output:
    tuple val(project), path(fastq), emit: combined_reads_projects

    script:
    project = "${proj_path}".split('/').last()

    """
    """

}

workflow {

    // Workflow for running the primer-trimming survey on sequencing reads

    // Probably a better way of doing this, keep for now because it works
    COMBINE_READS_PROJECTS(reads, projects) // Need to track both project id and paired-end reads

    // QC
    // fastp(COMBINE_READS_PROJECTS.out.combined_reads_projects, outdir)
    // seqprep_merge(fastp.out.cleaned_fastq, outdir)
    
    // Check for presence of primers
    STD_PRIMER_FLAG(COMBINE_READS_PROJECTS.out.combined_reads_projects, outdir) // Standard Library primers
    // general_primer_flag(seqprep_merge.out.merged_fastq, outdir) // Primers in general

    // Combining std and general primer outputs and parsing them to guide 
    // samples through automatic primer identification and trimming by cutadapt
    // comb_flags = general_primer_flag.out.general_primer_out
    // .join(std_primer_flag.out.std_primer_out)
    // trimming_conductor(comb_flags, outdir) // Generate a flags file with vals of 'none/std/auto' for both fwd and rev
    // parse_conductor(trimming_conductor.out.trimming_conductor_out, outdir) // Parse flags file into env variables
    
    // Join flags with merged fastq files
    // auto_trimming_input = parse_conductor.out.conductor_out
    // .join(seqprep_merge.out)

    // Run subworkflow for automatic primer trimming
    // Outputs either empty fasta file or fasta file containing predicted primers
    // AUTOMATIC_PRIMER_TRIMMING(
    //     auto_trimming_input,
    //     outdir
    //     )

    // Join auto flags to std flags and generated a concatenated fasta file containing primers to trim off
    // This can contain any valid combination of stranded std/auto primers 
    // concat_input = std_primer_flag.out.std_primer_out
    // .join(AUTOMATIC_PRIMER_TRIMMING.out.auto_primer_trimming_out)
    // concat_primers(concat_input, outdir)

    // Join concatenated primers to the fastp-cleaned paired reads files and run cutadapt on them
    // cutadapt_input = concat_primers.out.concat_primers_out
    // .join(fastp.out.cleaned_fastq)
    // cutadapt(cutadapt_input, outdir)

    // Just some logging for myself, will delete this eventually
    // final_out = fastp.out.cleaned_fastq
    // .join(cutadapt.out.cutadapt_out, remainder: true)
    // .map( { if (it[4] == null) { tuple(it[0], it[1], it[2], it[3]) } else { tuple(it[0], it[1], it[4], it[5]) }} )


}