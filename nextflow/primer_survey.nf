nextflow.enable.dsl=2

// nextflow run nextflow/primer_survey.nf --path /hps/nobackup/rdf/metagenomics/service-team/users/chrisata/asv_nf_testing/data --outdir /hps/nobackup/rdf/metagenomics/service-team/users/chrisata/asv_nf_testing
// nextflow run nextflow/primer_survey.nf -profile lsf --path /hps/nobackup/rdf/metagenomics/service-team/users/chrisata/asv_datasets/survey_datasets/final_survey_dataset_data/data --outdir /hps/nobackup/rdf/metagenomics/service-team/users/chrisata/asv_datasets/survey_datasets/final_survey_dataset_data/merged

include { fastp } from './modules/fastp.nf'
include { seqprep_merge } from './modules/seqprep_merge.nf'
include { std_primer_flag } from './modules/std_primer_flag.nf'
include { general_primer_flag } from './modules/general_primer_flag.nf'
include { trimming_conductor } from './modules/trimming_conductor.nf'
include { parse_conductor } from './modules/parse_conductor.nf'

include { automatic_primer_trimming } from './subworkflows/automatic_primer_trimming.nf'


params.path = null
params.outdir = null
params.cpu = null
params.filter = file("$projectDir/assets/NO_FILE")

reads = Channel.fromFilePairs( "${params.path}/DRP**/raw/*_{1,2}.fastq.gz" )
// reads = Channel.fromPath( "${params.outdir}/merged/**/*_MERGED.fastq.gz" )
projects = Channel.fromPath( "${params.path}/DRP*", type: 'dir' )
outdir = params.outdir
cpu = params.cpu

process COMBINE_READS_PROJECTS {

    label 'light'

    input:
    tuple val(sampleId), path(fastq)
    // path fastq

    val proj_path

    output:
    tuple val(sampleId), path(fastq), val(project), emit: combined_reads_projects
    // tuple path(fastq), val(project), emit: combined_reads_projects


    script:
    project = "${proj_path}".split('/').last()

    """
    """

}

process std_trimmer {

    label 'light'
    
    input:
    tuple val(project), val(fwd_flag), val(rev_flag)
    path std_primer_out
    path fastq_1
    path fastq_2

    output:
    stdout

    """
    if [[ ${fwd_flag} = "std" ]] | [[ ${rev_flag} = "std" ]]; then
        echo 'std trimming!'
    else
        echo 'not std trimming!'
    fi
    """

}

process cutadapt {

    label 'light'
    
    input:
    tuple val(project), path(std_primers), path(auto_primers)

    output:
    stdout

    // script:
    // def auto_filter = auto_primers.name != "NO_FILE" ? "$auto_primers" : ''

    """
    if [[ -s $auto_primers ]]; then
        echo $auto_primers
    else
        echo 'no auto primers!'
    fi
    """
}


workflow {

    COMBINE_READS_PROJECTS(reads, projects)
    fastp(COMBINE_READS_PROJECTS.out.combined_reads_projects, outdir)
    seqprep_merge(fastp.out.cleaned_fastq, outdir)
    std_primer_flag(seqprep_merge.out.merged_fastq, outdir)
    general_primer_flag(seqprep_merge.out.merged_fastq, outdir)

    // general_primer_flag(COMBINE_READS_PROJECTS.out.combined_reads_projects, outdir)

    comb_flags = general_primer_flag.out.general_primer_out
    .join(std_primer_flag.out.std_primer_out)
    
    trimming_conductor(comb_flags, outdir)
    parse_conductor(trimming_conductor.out.trimming_conductor_out, outdir)
    
    // std_trimmer(
    //     parse_conductor.out.conductor_out,
    //     std_primer_flag.out.std_primer_out.map{ it[1] }, 
    //     fastp.out.cleaned_fastq.map{ it[1] }, 
    //     fastp.out.cleaned_fastq.map{ it[2] }
    //     ).collect()
    //     .view()
    
    auto_trimming_input = parse_conductor.out.conductor_out.join(seqprep_merge.out)
    automatic_primer_trimming(
        auto_trimming_input,
        outdir
        )

    cutadapt_input = std_primer_flag.out.std_primer_out
    .join(automatic_primer_trimming.out.auto_primer_trimming_out, remainder: true)

    null_filter = cutadapt_input.map( {if (it[2] == null){ tuple(it[0], params.filter) } else { tuple(it[0], it[2]) } } )

    final_cutadapt_input = std_primer_flag.out.std_primer_out
    .join(null_filter)

    cutadapt(final_cutadapt_input).view()

}