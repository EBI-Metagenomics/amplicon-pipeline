
include { ASSESS_MCP_CONS } from '../modules/assess_mcp_cons.nf'
include { FIND_MCP_INF_POINTS } from '../modules/find_mcp_inf_points.nf'
include { ASSESS_MCP_INF_POINTS } from '../modules/assess_mcp_inf_points.nf'

workflow AUTOMATIC_PRIMER_PREDICTION {

    // Subworkflow that performs automatic primer trimming on sequencing reads
    // Takes flags for the forward and reverse strands to make sure primers are  
    // predicted only for the required strands
    // In theory this should work on single-end reads as 'forward' and 'reverse'
    // more accurately means 5' or 3' ends of a merged paired-end file
    // Outputs a fasta file that is either empty or that contains the predicted primer sequences (5' to 3')
    
    take:
        auto_trimming_input
        outdir
    main:
        // Use Most Common Prefix (MCP) method to generate curves of base conservation
        ASSESS_MCP_CONS(
            auto_trimming_input.map{ it[0] }, // project ID
            auto_trimming_input.map{ it[1] }, // Sample ID
            auto_trimming_input.map{ it[2] }, // fwd_flag
            auto_trimming_input.map{ it[3] }, // rev_flag
            auto_trimming_input.map{ it[4] }, // fastq
            outdir
        )

        // Find inflection points in conservation curves
        FIND_MCP_INF_POINTS(
            ASSESS_MCP_CONS.out.mcp_cons_out,
            outdir
        ) 

        // Join fastq channel and the inf_points channel
        assess_inf_input = FIND_MCP_INF_POINTS.out.inf_points_out 
                           .join(auto_trimming_input.map{ it[0, 1, 4] }, by: [0, 1])
        
        // Select inflection points most likely to be primer cutoff points
        ASSESS_MCP_INF_POINTS(
            assess_inf_input,
            outdir
        )

   emit:
        // Had to make even std runs go through this workflow (albeit just outputting empty files for every process) due to a join later being slow otherwise. 
        // wonder if there's a better way?
        auto_primer_trimming_out = ASSESS_MCP_INF_POINTS.out

}