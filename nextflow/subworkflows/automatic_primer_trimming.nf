
include { assess_mcp_cons } from '../modules/assess_mcp_cons.nf'
include { find_mcp_inf_points } from '../modules/find_mcp_inf_points.nf'
include { assess_mcp_inf_points } from '../modules/assess_mcp_inf_points.nf'

workflow AUTOMATIC_PRIMER_TRIMMING {

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
        assess_mcp_cons(
            auto_trimming_input.map{ it[0] }, // project ID
            auto_trimming_input.map{ it[1] }, // fwd_flag
            auto_trimming_input.map{ it[2] }, // rev_flag
            auto_trimming_input.map{ it[3] }, // fastq
            outdir
            )

        find_mcp_inf_points(assess_mcp_cons.out.mcp_cons_out, outdir) // Find inflection points in conservation curves
        assess_inf_input = find_mcp_inf_points.out.inf_points_out // Join fastq channel and the inf_points channel
        .join(auto_trimming_input.map{ it[0, 3] })
        assess_mcp_inf_points(assess_inf_input, outdir) // Select inflection points most likely to be primer cutoff points

   emit:
        // Had to make even std runs go through this workflow (albeit just outputting empty files for every process) due to a join later being slow otherwise. 
        // wonder if there's a better way?
        auto_primer_trimming_out = assess_mcp_inf_points.out

}