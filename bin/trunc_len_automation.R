
# Heavily inspired by what @alexandriai168  wrote on this GitHub issue:
# https://github.com/benjjneb/dada2/issues/1729

trunc_len_automation = function(fastq_path){
    
    n = 500000
    reads_remainder_threshold = 0.5
    final_where_to_cut = 0

    srqa = qa(fastq_path, n=n)
    df = srqa[["perCycle"]]$quality
    
    means = rowsum(df$Score*df$Count, df$Cycle)/rowsum(df$Count, df$Cycle)
    total_read_count = srqa[["readCounts"]]$read
    counts_per_cycle = rowsum(df$Count, df$Cycle)

    q_score_list = 20:30

    for (q in q_score_list){
        less_than_q = min(which(means<q)) - 1
        reads_remainder = counts_per_cycle[less_than_q]/total_read_count
        
        if (less_than_q != Inf){
            if (reads_remainder >= reads_remainder_threshold){
                final_where_to_cut = less_than_q
                break
            }
        }
    }

    return(final_where_to_cut)
}