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


# Have to use `box` instead of `library` and `source` so that custom scripts can be loaded when executed by Nextflow
box::use(tidyverse[...])
box::use(data.table[...])
box::use(dada2[...])

# Custom function for tracking reads to their DADA2-generated ASVs (bin/read_asv_tracking.R) 
box::use(./read_asv_tracking[...])
# Custom function for automatic truncation of reads based on quality scores (bin/trunc_len_automation.R)
box::use(./trunc_len_automation[...])

args = commandArgs(trailingOnly=TRUE) # Expects at most three arguments, a prefix, and one fastq for each strand (F and R)
                                      # If it's a single-end run, then the third argument should not be used
prefix = args[1] # Prefix
path_f = args[2] # Forward fastq
path_r = args[3] # Reverse fastq

# different tax ranks for silva/pr2
silva_tax_vec = c("Superkingdom", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
pr2_tax_vec = c("Domain", "Supergroup", "Division", "Subdivision", "Class", "Order", "Family", "Genus", "Species")

# Identify truncLen parameter for filterAndTrim function
final_where_to_cut_f = trunc_len_automation(path_f)
if (!is.na(path_r)){
  final_where_to_cut_r = trunc_len_automation(path_r)
}

# Do some quality filtering
filt_f =  paste0("./", prefix, "_1", "_filt.fastq.gz")
tryCatch(
  {
    if (!is.na(path_r)){
      filt_r =  paste0("./", prefix, "_2", "_filt.fastq.gz")
      print(paste0("The forward strand truncation point is: ", final_where_to_cut_f))
      print(paste0("The reverse strand truncation point is: ", final_where_to_cut_r))
      out = filterAndTrim(path_f, filt_f, path_r, filt_r, rm.phix=TRUE, maxEE=c(2,5), truncQ=2, truncLen=c(final_where_to_cut_f,final_where_to_cut_r), compress=TRUE, multithread=TRUE)
    } else{
      print(paste0("The forward strand truncation point is: ", final_where_to_cut_f))
      out = filterAndTrim(path_f, filt_f, rm.phix=TRUE, maxEE=2, truncQ=2, truncLen=final_where_to_cut_f, compress=TRUE, multithread=TRUE)
    }
  }, error = function(msg){
    message(paste("Caught an error at the `filterAndTrim` stage:\n", msg))
    quit()
  }
)

tryCatch(
  {
    # Learn error model
    err_f = learnErrors(filt_f, multithread=TRUE)
    if (!is.na(path_r)){
      err_r = learnErrors(filt_r, multithread=TRUE)
    }
  }, error = function(msg){
    message(paste("Caught an error at the `learnErrors` stage:\n", msg))
    quit()
  }
)

tryCatch(
  {
    # Dereplicate sequences
    drp_f = derepFastq(filt_f)
    if (!is.na(path_r)){
      drp_r = derepFastq(filt_r)
    }
  }, error = function(msg){
    message(paste("Caught an error at the `derepFastq` stage:\n", msg))
    quit()
  }
)

tryCatch(
  {
    # Generate stranded ASVs
    dada_f = dada(drp_f, err=err_f, multithread=TRUE)
    if (!is.na(path_r)){
      dada_r = dada(drp_r, err=err_r, multithread=TRUE)
    }
  }, error = function(msg){
    message(paste("Caught an error at the `dada` stage:\n", msg))
    quit()
  }
)

tryCatch(
  {
    # Merge stranded ASVs
    if (!is.na(path_r)){
      merged = mergePairs(dada_f, drp_f, dada_r, drp_r, verbose=TRUE)
    } else{
      merged = dada_f
    }
  }, error = function(msg){
    message(paste("Caught an error at the `mergePairs` stage:\n", msg))
    quit()
  }
)

if (length(merged$sequence) == 0){
  message("Caught an error - No ASVs - stopping script early.")
  quit()
}

tryCatch(
  {
    # Make ASV count table
    seqtab = makeSequenceTable(merged)
  }, error = function(msg){
    message(paste("Caught an error at the `makeSequenceTable` stage:\n", msg))
    quit()
  }
)

tryCatch(
  {
    # Remove chimeras
    seqtab.nochim = removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
  }, error = function(msg){
    message(paste("Caught an error at the `removeBimeraDenovo` stage:\n", msg))
    quit()
  }
)

chimera_ids = which(colnames(seqtab) %in% colnames(seqtab.nochim) == FALSE)

tryCatch(
  {
    # Track reads to their ASVs
    if (!is.na(path_r)){
      final_f_map = read_asv_tracking(dada_f, drp_f, merged, "forward", chimera_ids)
      final_r_map = read_asv_tracking(dada_r, drp_r, merged, "reverse", chimera_ids)
    } else{
      final_f_map = read_asv_tracking(dada_f, drp_f, merged, "single", chimera_ids)
    }
  }, error = function(msg){
    message(paste("Caught an error at the `read_asv_tracking` stage:\n", msg))
    quit()
  }
)

final_f_output = list()
if (!is.na(path_r)){
  final_r_output = list()
}

unmatched_asvs = character()

for (i in 1:length(final_f_map)){
  f_map_list = final_f_map[[i]]
  if (!is.na(path_r)){
    r_map_list = final_r_map[[i]]
    overlap = intersect(f_map_list, r_map_list)
    
    if (length(overlap) == 0){
      unmatched_asvs = append(unmatched_asvs, i)
    } else if (overlap == 0){
      unmatched_asvs = append(unmatched_asvs, i)
    }
  }

  final_f_output[[i]] = f_map_list[1]
  if (!is.na(path_r)){
    final_r_output[[i]] = f_map_list[1]
  }

}
# The extremely vast majority of forwards+reverse pairs should be assigned the same ASV. This checks it
traced_remainder = length(final_f_map) - length(unmatched_asvs)
total_dada2_reads = sum(seqtab.nochim)
final_matched_perc = traced_remainder / total_dada2_reads

# Save map to file, with each line representing an ASV's read
if (!is.na(path_r)){
  fwrite(final_f_output, file = paste0("./", prefix, "_1_map.txt"), sep="\n")
  fwrite(final_r_output, file = paste0("./", prefix, "_2_map.txt"), sep="\n")
} else{
  fwrite(final_f_output, file = paste0("./", prefix, "_map.txt"), sep="\n")
}

# Save ASV count table
write.table(seqtab.nochim, file = paste0("./", prefix, "_asv_counts.tsv"), sep = "\t", row.names=FALSE)

# Write proportion of chimeric reads into a file
seqtab_read_count = sum(seqtab)
seqtab.nochim_read_count = sum(seqtab.nochim)
proportion_chimeric = 1 - (seqtab.nochim_read_count / seqtab_read_count)

# Get count of unique ASVs left after all types of filtering
asvs_left = sort(as.numeric(unique(unlist(final_f_output))))
asvs_left = asvs_left[2:length(asvs_left)]

# Save ASV sequences to FASTA file
seqtab.length = length(seqtab)
num_list = as.character(1:seqtab.length)
id_list = paste("seq", asvs_left, sep="_")
unqs = getUniques(seqtab)[asvs_left]
uniquesToFasta(unqs, paste0("./", prefix, "_asvs.fasta"), id_list)

output_report_df <- data.frame(
  names = c("initial_number_of_reads", "proportion_matched", "proportion_chimeric", "final_number_of_reads"),
  values = c(length(final_f_map), final_matched_perc, proportion_chimeric, total_dada2_reads)
)
write.table(output_report_df, file = paste0("./", prefix, "_dada2_stats.tsv"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
