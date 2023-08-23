
library(tidyverse)
library(dada2)
library(data.table)

# Load function for tracking reads to their DADA2-generated ASVs
source("/hps/software/users/rdf/metagenomics/service-team/users/chrisata/asv_gen/bin/read_asv_tracking.R")

args = commandArgs(trailingOnly=TRUE) # Expects thre arguments, one fastq for each strand (F and R), and a prefix

path_f = args[1] # Forward fastq
path_r = args[2] # Reverse fastq
prefix = args[3] # Prefix

# Learn error model
err_f = learnErrors(path_f, multithread=TRUE)
err_r = learnErrors(path_r, multithread=TRUE)

# Dereplicate sequences
drp_f = derepFastq(path_f)
drp_r = derepFastq(path_r)

# Generate stranded ASVs
dada_f = dada(drp_f, err=err_f, multithread=TRUE)
dada_r = dada(drp_r, err=err_r, multithread=TRUE)

# Merge stranded ASVs
merged = mergePairs(dada_f, drp_f, dada_r, drp_r, verbose=TRUE)

# Make ASV count table
seqtab = makeSequenceTable(merged)

# Remove chimeras
seqtab.nochim = removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

# Assign taxonomy (using SILVA)
taxa = assignTaxonomy(seqtab.nochim, "./data/silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)

chimera_ids = which(colnames(seqtab) %in% colnames(seqtab.nochim) == FALSE)

# Track reads to their ASVs
final_f_map = read_asv_tracking(dada_f, drp_f, merged, "forward", chimera_ids)
final_r_map = read_asv_tracking(dada_r, drp_r, merged, "reverse", chimera_ids)

final_f_output = list()
final_r_output = list()

unmatched_asvs = character()

for (i in 1:length(final_f_map)){
  f_map_list = final_f_map[[i]]
  r_map_list = final_r_map[[i]]
  
  overlap = intersect(f_map_list, r_map_list)  
  
  if (length(overlap) == 0){
    unmatched_asvs = append(unmatched_asvs, i)
  }
  else if (overlap == 0){
    unmatched_asvs = append(unmatched_asvs, i)
  }
  
  final_f_output[[i]] = paste(f_map_list, collapse = ",")
  final_r_output[[i]] = paste(r_map_list, collapse = ",")
  
}


# The extremely vast majority of forwards+reverse pairs should be assigned the same ASV. This checks it
traced_remainder = length(final_f_map) - length(unmatched_asvs)
total_dada2_reads = sum(seqtab.nochim)
final_matched_perc = traced_remainder / total_dada2_reads
write(final_matched_perc, paste0("./", prefix, "_proportion_matched.txt"))

# Save map to file, with each line representing an ASV's read 
fwrite(final_f_output, file = paste0("./", prefix, "_1_map.txt"), sep="\n")
fwrite(final_r_output, file = paste0("./", prefix, "_2_map.txt"), sep="\n")

# Save taxa annotation to file
taxa = cbind(rownames(taxa), taxa)
colnames(taxa) = c("ASV", colnames(taxa)[2:7])
write.table(taxa, file = paste0("./", prefix, "_taxa.tsv"), sep = "\t", row.names=FALSE)

# Write proportion of chimeric reads into a file
seqtab_read_count = sum(seqtab)
seqtab.nochim_read_count = sum(seqtab.nochim)
proportion_chimeric = 1 - (seqtab.nochim_read_count / seqtab_read_count)
write(proportion_chimeric, paste0("./", prefix, "_proportion_chimeric.txt"))
