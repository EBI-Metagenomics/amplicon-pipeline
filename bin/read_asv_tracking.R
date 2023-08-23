
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

  mergers_map = list()
  
  for (i in 1:length(dada_obj$sequence)) {
    if ((i %in% merged.strand)){
      matches = which(merged.strand == i)
      add_lst = list()

      for (match in matches){

        if (!match %in% chimera_ids){
          add_lst = append(add_lst, match)
        }
        
      }
      
      if (length(add_lst) != 0){
        mergers_map[[i]] = add_lst
      }
      else{
        mergers_map[[i]] = list(0)
      }
            
    }
    else{
      mergers_map[[i]] = list(0)
    }

  }
  
  # Create map from stranded ASV to dereplicated sequences
  derep_map = list()
  counter = 0
  for (old_asv in dada_obj$map) {
    counter = counter + 1

    if (is.na(old_asv)){
      derep_map[[counter]] = list(0)
      next
    }
    
    matches = mergers_map[[old_asv]]
    add_lst = list()

    for (match in matches){
      if (match != 0){
        add_lst = append(add_lst, match)
      }
      
    }

    if (length(add_lst) != 0){
      derep_map[[counter]] = add_lst
    }
    else{
      derep_map[[counter]] = list(0)
    }

  }

  # Create map from dereplicated sequences to all sequences
  # final_map = rep("0,", length(drp_obj$map))
  final_map = list()
  
  counter = 0
  for (old_derep in drp_obj$map) {
    counter = counter + 1

    if (is.na(old_derep)){
      final_map[[counter]] = c(0)
      next
    }

    matches = derep_map[[old_derep]]
    add_lst = character()
    
    for (match in matches){
      if (match != 0){
        add_lst = append(add_lst, match)
      }
      
    }

    if (length(add_lst) != 0){
      final_map[[counter]] = as.numeric(add_lst)
    }
    else{
      final_map[[counter]] = c(0)
    }  
  
  }

  return(final_map)
}