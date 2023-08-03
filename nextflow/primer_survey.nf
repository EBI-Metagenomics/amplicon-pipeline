nextflow.enable.dsl=2

// nextflow run nextflow/primer_survey.nf --path /hps/nobackup/rdf/metagenomics/service-team/users/chrisata/asv_nf_testing/data --outdir /hps/nobackup/rdf/metagenomics/service-team/users/chrisata/asv_nf_testing
// nextflow run nextflow/primer_survey.nf -profile lsf --path /hps/nobackup/rdf/metagenomics/service-team/users/chrisata/asv_datasets/survey_datasets/final_survey_dataset_data/data --outdir /hps/nobackup/rdf/metagenomics/service-team/users/chrisata/asv_datasets/survey_datasets/final_survey_dataset_data/merged

include { fastp } from './modules/fastp.nf'
include { seqprep_merge } from './modules/seqprep_merge.nf'
include { std_primer_flag } from './modules/std_primer_flag.nf'
include { general_primer_flag } from './modules/general_primer_flag.nf'
include { trimming_conductor } from './modules/trimming_conductor.nf'
include { parse_conductor } from './modules/parse_conductor.nf'
include { cutadapt } from './modules/cutadapt.nf'
include { concat_primers } from './modules/concat_primers.nf'

include { automatic_primer_trimming } from './subworkflows/automatic_primer_trimming.nf'


params.path = null
params.outdir = null
params.cpu = null
params.filter = file("$projectDir/assets/NO_FILE")

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
    tuple val(project), val(sampleId), path(fastq), emit: combined_reads_projects

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

    comb_flags = general_primer_flag.out.general_primer_out
    .join(std_primer_flag.out.std_primer_out)
    

    trimming_conductor(comb_flags, outdir)
    parse_conductor(trimming_conductor.out.trimming_conductor_out, outdir)
    
    auto_trimming_input = parse_conductor.out.conductor_out.join(seqprep_merge.out)
    automatic_primer_trimming(
        auto_trimming_input,
        outdir
        )

    concat_input = std_primer_flag.out.std_primer_out
    .join(automatic_primer_trimming.out.auto_primer_trimming_out)

    concat_primers(concat_input, outdir)

    cutadapt_input = concat_primers.out.concat_primers_out
    .join(fastp.out.cleaned_fastq)
    
    cutadapt(cutadapt_input, outdir)

    final_out = fastp.out.cleaned_fastq
    .join(cutadapt.out.cutadapt_out, remainder: true)
    .map( { if (it[4] == null) { tuple(it[0], it[1], it[2], it[3]) } else { tuple(it[0], it[1], it[4], it[5]) }} )


}