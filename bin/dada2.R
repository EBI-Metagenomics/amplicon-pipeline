
library(tidyverse)
library(dada2)
library(data.table)
library(ShortRead)

# Load function for tracking reads to their DADA2-generated ASVs
source("/hps/software/users/rdf/metagenomics/service-team/users/chrisata/asv_gen/bin/read_asv_tracking.R")
source("/hps/software/users/rdf/metagenomics/service-team/users/chrisata/asv_gen/bin/trunc_len_automation.R")

args = commandArgs(trailingOnly=TRUE) # Expects thre arguments, one fastq for each strand (F and R), and a prefix

prefix = args[1] # Prefix
ref_label = args[2] # Reference DB name
ref_db = args[3] # Reference DB
path_f = args[4] # Forward fastq
path_r = args[5] # Reverse fastq

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
if (!is.na(path_r)){
  filt_r =  paste0("./", prefix, "_2", "_filt.fastq.gz")
  print(final_where_to_cut_f)
  print(final_where_to_cut_r)
  out = filterAndTrim(path_f, filt_f, path_r, filt_r, rm.phix=TRUE, maxEE=c(1,1), truncQ=2, truncLen=c(final_where_to_cut_f,final_where_to_cut_r), compress=TRUE, multithread=TRUE)
} else{
  print(final_where_to_cut_f)
  out = filterAndTrim(path_f, filt_f, rm.phix=TRUE, maxEE=1, truncQ=2, truncLen=final_where_to_cut_f, compress=TRUE, multithread=TRUE)
}

# Learn error model
err_f = learnErrors(filt_f, multithread=TRUE)
if (!is.na(path_r)){
  err_r = learnErrors(filt_r, multithread=TRUE)
}

# Dereplicate sequences
drp_f = derepFastq(filt_f)
if (!is.na(path_r)){
  drp_r = derepFastq(filt_r)
}

# Generate stranded ASVs
dada_f = dada(drp_f, err=err_f, multithread=TRUE)
if (!is.na(path_r)){
  dada_r = dada(drp_r, err=err_r, multithread=TRUE)
}

# Merge stranded ASVs

if (!is.na(path_r)){
  merged = mergePairs(dada_f, drp_f, dada_r, drp_r, verbose=TRUE)
} else{
  merged = dada_f
}

if (length(merged$sequence) == 0){
  print("No ASVs - stopping script early.")
} else{
    # Make ASV count table
  seqtab = makeSequenceTable(merged)

  # Remove chimeras
  seqtab.nochim = removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

  # Assign taxonomy (using SILVA)
  if (ref_label == "DADA2-SILVA"){
    taxa = assignTaxonomy(seqtab.nochim, ref_db, multithread=TRUE, taxLevels=silva_tax_vec)
  } else if (ref_label == "DADA2-PR2"){
    taxa = assignTaxonomy(seqtab.nochim, ref_db, multithread=TRUE, taxLevels=pr2_tax_vec)
  }

  chimera_ids = which(colnames(seqtab) %in% colnames(seqtab.nochim) == FALSE)

  # Track reads to their ASVs
  if (!is.na(path_r)){
    final_f_map = read_asv_tracking(dada_f, drp_f, merged, "forward", chimera_ids)
    final_r_map = read_asv_tracking(dada_r, drp_r, merged, "reverse", chimera_ids)
  } else{
    final_f_map = read_asv_tracking(dada_f, drp_f, merged, "single", chimera_ids)
  }

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
    
    final_f_output[[i]] = paste(f_map_list, collapse = ",")
    if (!is.na(path_r)){
      final_r_output[[i]] = paste(r_map_list, collapse = ",")
    }
    
  }

  # The extremely vast majority of forwards+reverse pairs should be assigned the same ASV. This checks it
  traced_remainder = length(final_f_map) - length(unmatched_asvs)
  total_dada2_reads = sum(seqtab.nochim)
  final_matched_perc = traced_remainder / total_dada2_reads
  write(final_matched_perc, paste0("./", prefix, "_proportion_matched.txt"))

  # Save map to file, with each line representing an ASV's read
  if (!is.na(path_r)){
    fwrite(final_f_output, file = paste0("./", prefix, "_1_map.txt"), sep="\n")
    fwrite(final_r_output, file = paste0("./", prefix, "_2_map.txt"), sep="\n")
  } else{
    fwrite(final_f_output, file = paste0("./", prefix, "_map.txt"), sep="\n")
  }
  # Save taxa annotation to file
  taxa = cbind(rownames(taxa), taxa)
  if (ref_label == "DADA2-SILVA"){
    colnames(taxa) = c("ASV", colnames(taxa)[2:9])
  } else if (ref_label == "DADA2-PR2"){
    colnames(taxa) = c("ASV", colnames(taxa)[2:10])
  }
  write.table(taxa, file = paste0("./", prefix, "_", ref_label, "_taxa.tsv"), sep = "\t", row.names=FALSE)

  # Save ASV count table
  write.table(seqtab.nochim, file = paste0("./", prefix, "_asv_counts.tsv"), sep = "\t", row.names=FALSE)


  # Write proportion of chimeric reads into a file
  seqtab_read_count = sum(seqtab)
  seqtab.nochim_read_count = sum(seqtab.nochim)
  proportion_chimeric = 1 - (seqtab.nochim_read_count / seqtab_read_count)
  write(proportion_chimeric, paste0("./", prefix, "_proportion_chimeric.txt"))


}
