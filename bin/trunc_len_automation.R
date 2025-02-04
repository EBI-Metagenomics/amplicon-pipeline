#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-

# Copyright 2024 EMBL - European Bioinformatics Institute
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

#' @export
trunc_len_automation = function(fastq_path){
  # Heavily inspired by what @alexandriai168 wrote on this GitHub issue:
  # https://github.com/benjjneb/dada2/issues/1729

  box::use(ShortRead[...])
  box::use(utils[tail])

  n = 500000
  reads_remainder_threshold = 0.5
  final_where_to_cut = 0
  
  srqa = qa(fastq_path, n=n)
  df = srqa[["perCycle"]]$quality
  
  total_read_count = srqa[["readCounts"]]$read
  counts_per_cycle = rowsum(df$Count, df$Cycle)
  means = tail(rowsum(df$Score*df$Count, df$Cycle)/rowsum(df$Count, df$Cycle), length(counts_per_cycle) - 99)
  filtered_counts_per_cycle = tail(counts_per_cycle, length(counts_per_cycle) - 99)
  
  q_score_list = 20:30
  
  for (q in q_score_list){
    less_than_q = min(which(means<q)) - 1
    
    if (less_than_q > 0 & less_than_q < length(filtered_counts_per_cycle)){

      reads_remainder = filtered_counts_per_cycle[less_than_q]/total_read_count
      if (less_than_q != Inf){
        if (reads_remainder >= reads_remainder_threshold){
          final_where_to_cut = less_than_q
          break
        }
      }
    }
  }

  if (final_where_to_cut == 0){
    return(0)
  }else {
    final_where_to_cut = as.integer(rownames(filtered_counts_per_cycle)[final_where_to_cut])
    return(final_where_to_cut)
  }
}
