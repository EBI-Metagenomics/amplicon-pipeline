#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { validateParameters; paramsHelp; paramsSummaryLog; fromSamplesheet } from 'plugin/nf-validation'

// Print help message, supply typical command line usage for the pipeline
if (params.help) {
   log.info paramsHelp("nextflow run main.nf --input samplesheet.csv --outdir ./outdir")
   exit 0
}
// Validate input parameters
validateParameters()

// Print summary of supplied parameters
log.info paramsSummaryLog(workflow)

ch_input = Channel.fromSamplesheet( "input" )

include { AMPLICON_PIPELINE_V6 } from './workflows/pipeline.nf'

workflow {
    AMPLICON_PIPELINE_V6()
}