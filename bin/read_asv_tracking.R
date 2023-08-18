
read_asv_tracking = function(dada_obj, drp_obj, merged, strand, chimera_ids) {
  
  # Select requested strand ASV map
  if (strand == "forward"){
    merged.strand = merged$forward
  }
  else if (strand == "reverse") {
    merged.strand = merged$reverse
  }
  
  if (length(chimera_ids) != 0){
    merged.strand = merged.strand[-chimera_ids]
  }
  
  # Create map from merged ASV to stranded ASV
  mergers_map = numeric(length(dada_obj$sequence))
  for (i in 1:length(dada_obj$sequence)) {
    
    if ((i %in% merged.strand)){
      matches = which(merged.strand == i)
      for (match in matches){
        if (!match %in% chimera_ids){
          mergers_map[i] = match
        }
        else{
          mergers_map[i] = NA
        }
      }
    }
    else {
      mergers_map[i] = NA
    }
  }

  
  # Create map from stranded ASV to dereplicated sequences
  derep_map = numeric(length(dada_obj$map))
  counter = 0
  for (i in dada_obj$map) {
    counter = counter + 1
    if (!is.na(mergers_map[i])){
      derep_map[counter] = mergers_map[i]
    }
    else {
      derep_map[counter] = NA
    }
  }
  
  # Create map from dereplicated sequences to all sequences
  final_map = numeric(length(drp_obj$map))
  counter = 0
  for (i in drp_obj$map) {
    counter = counter + 1
    if (!is.na(derep_map[i])){
      final_map[counter] = derep_map[i]
    }
    else {
      final_map[counter] = 0
    }
  }
  
  return(final_map)

}