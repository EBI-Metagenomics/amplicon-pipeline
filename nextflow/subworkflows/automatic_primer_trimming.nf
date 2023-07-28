
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
        // .collect()
        // .view()

        find_mcp_inf_points(assess_mcp_cons.out.mcp_cons_out, outdir)
        // .collect()
        // .view()

        assess_inf_input = find_mcp_inf_points.out.inf_points_out.join(auto_trimming_input.map{ it[0, 3] })

        assess_mcp_inf_points(assess_inf_input, outdir)
        // .collect()
        // .view()
   emit:
     auto_primer_trimming_out = assess_mcp_inf_points.out

}