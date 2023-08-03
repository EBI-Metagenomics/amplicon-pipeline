
include { assess_mcp_cons } from '../modules/assess_mcp_cons.nf'
include { find_mcp_inf_points } from '../modules/find_mcp_inf_points.nf'
include { assess_mcp_inf_points } from '../modules/assess_mcp_inf_points.nf'

workflow automatic_primer_trimming {
    
    // [project, fwd_flag, rev_flag, fastq, outdir]
    take:
        auto_trimming_input
        outdir
    main:
        assess_mcp_cons(
            auto_trimming_input.map{ it[0] },
            auto_trimming_input.map{ it[1] },
            auto_trimming_input.map{ it[2] },
            auto_trimming_input.map{ it[3] },
            outdir
            )

        find_mcp_inf_points(assess_mcp_cons.out.mcp_cons_out, outdir)
        assess_inf_input = find_mcp_inf_points.out.inf_points_out.join(auto_trimming_input.map{ it[0, 3] })
        assess_mcp_inf_points(assess_inf_input, outdir)

   emit:
//    Had to make even std runs go through this workflow (albeit just outputting empty files for every process) due to a join later being slow otherwise. 
//    Wonder if there's a better way?
     auto_primer_trimming_out = assess_mcp_inf_points.out

}