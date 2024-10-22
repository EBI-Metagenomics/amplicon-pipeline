
include { ASSESSMCPPROPORTIONS   } from '../../modules/ebi-metagenomics/assessmcpproportions/main'
include { FIND_MCP_INF_POINTS   } from '../../modules/local/find_mcp_inf_points/main.nf'
include { ASSESS_MCP_INF_POINTS } from '../../modules/local/assess_mcp_inf_points/main.nf'

workflow AUTOMATIC_PRIMER_PREDICTION {

    // Subworkflow that performs automatic primer trimming on sequencing reads
    // Takes flags for the forward and reverse strands to make sure primers are  
    // predicted only for the required strands
    // In theory this should work on single-end reads as 'forward' and 'reverse'
    // more accurately means 5' or 3' ends of a merged paired-end file
    // Outputs a fasta file that is either empty or that contains the predicted primer sequences (5' to 3')
    
    take:
        auto_trimming_input
    main:

        ch_versions = Channel.empty()
        // Use Most Common Prefix (MCP) method to generate curves of base conservation
        ASSESSMCPPROPORTIONS(
            auto_trimming_input,
            false
        )
        ch_versions = ch_versions.mix(ASSESSMCPPROPORTIONS.out.versions.first())

        // Find inflection points in conservation curves
        FIND_MCP_INF_POINTS(
            ASSESSMCPPROPORTIONS.out.tsv
        )
        ch_versions = ch_versions.mix(FIND_MCP_INF_POINTS.out.versions.first())

        extracted_reads = auto_trimming_input.map{ meta, fwd_flag, rev_flag, extracted_reads ->
                                                    [ meta, extracted_reads ]
                                                 }
        // Join fastq channel and the inf_points channel
        assess_inf_input = FIND_MCP_INF_POINTS.out.inf_points_out
                           .join(extracted_reads, by: [0])

        // Select inflection points most likely to be primer cutoff points
        ASSESS_MCP_INF_POINTS(
            assess_inf_input
        )
        ch_versions = ch_versions.mix(ASSESS_MCP_INF_POINTS.out.versions.first())

        ASSESS_MCP_INF_POINTS.out.auto_primer_out

   emit:
        auto_primer_trimming_out = ASSESS_MCP_INF_POINTS.out.auto_primer_out
        versions = ch_versions

}