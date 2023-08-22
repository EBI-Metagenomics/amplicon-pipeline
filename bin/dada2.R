
library(tidyverse)
library(dada2)
library(data.table)

# Load function for tracking reads to their DADA2-generated ASVs
source("./lib/read_asv_tracking.R")

args = commandArgs(trailingOnly=TRUE) # Expects thre arguments, one fastq for each strand (F and R), and a prefix

# path_f = args[1] # Forward fastq
# path_r = args[2] # Reverse fastq

path_f = "./MiSeq_SOP/filtered/F3D2_F_filt.fastq.gz"
path_r = "./MiSeq_SOP/filtered/F3D2_R_filt.fastq.gz"
prefix = "F3D2"

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
taxa = assignTaxonomy(seqtab.nochim, "./silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)

chimera_ids = which(colnames(seqtab) %in% colnames(seqtab.nochim) == FALSE)

# Track reads to their ASVs
final_f_map = read_asv_tracking(dada_f, drp_f, merged, "forward", chimera_ids)
final_r_map = read_asv_tracking(dada_r, drp_r, merged, "reverse", chimera_ids)

# The extremely vast majority of forwards+reverse pairs should be assigned the same ASV. This checks it
unmatched_asvs = which((final_r_map == final_f_map) == FALSE)
unmatched_asvs = unique(c(which(final_f_map == 0), which(final_r_map == 0), unmatched_asvs))
unmatched_perc = length(unmatched_asvs) / length(final_f_map) * 100

traced_remainder = length(final_f_map) - length(unmatched_asvs)
total_dada2_reads = sum(seqtab.nochim)
final_matched_perc = traced_remainder / total_dada2_reads * 100

# Save map to file, with each line representing an ASV's read 
fwrite(list(final_f_map), file = paste0("./", prefix, "_1_map.txt"))
fwrite(list(final_r_map), file = paste0("./", prefix, "_2_map.txt"))
taxa = cbind(rownames(taxa), taxa)
colnames(taxa) = c("ASV", colnames(taxa)[2:7])

write.table(taxa, file = paste0("./", prefix, "_taxa.tsv"), sep = "\t", row.names=FALSE)

# Write proportion of chimeric reads into a file
seqtab_read_count = sum(seqtab)
seqtab.nochim_read_count = sum(seqtab.nochim)
proportion_chimeric = 1 - (seqtab.nochim_read_count / seqtab_read_count)
write(proportion_chimeric, "./proportion_chimeric.txt")


