#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { AMPLICON_PIPELINE_V6 } from './workflows/pipeline.nf'

workflow {
    AMPLICON_PIPELINE_V6()
}