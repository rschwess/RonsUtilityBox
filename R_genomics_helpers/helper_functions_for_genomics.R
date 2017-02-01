# Various helper functions for genomics and related analysis
# Ron Schwessinger

# MAKE A RANDOM SEQUENCE WITH SEQ OR PFM MATCHES ADDED ------------------------
makeRandomSequence <- function(l, add.sequence = FALSE, add.pfm.match = FALSE, add.reverse=FALSE){
  # Input: 
  #   l = target sequence length
  #   add.sequence = default FALSE add a provided sequence (charstring) or list of charstrings, randomly to the sequence
  #   add.pfm.match = default FALSE, sample and add a sampled from the provided Pfm matrix frequencies, can be a list of pfms
  #   add.reverse = default FALSE, if adding a sequence or Pfm match randomly decide if to add forward or reverse complement
  # Returns: char string of random sequence of length l
  
  require(Biostrings)
  
  s <- paste(sample(c("A","C","G","T"), l, replace=TRUE), collapse="", sep="")

  # add Sequence 
  if("list" %in% is(add.sequence) || add.sequence != FALSE){
    
    # if(!"list" %in% is(add.sequence)){  # make into temp list if single sequence
    #   add.sequence <- list(add.sequence)
    # }
    # 
    for(i in c(1:length(add.sequence))){  # add a sequence for every entrie in the list
      add.sequence[[i]] <- as.character(add.sequence[[i]])
      offset <- sample(c(1:(l-nchar(add.sequence[[i]]))), 1) # sample offset where to add
      if(!add.reverse){
        # add straight forward
        substr(s,offset,offset+nchar(add.sequence[[i]])) <- add.sequence[[i]]
      }else if(add.reverse){
        #randomly decide if to add forward or reverse complement
        j <- sample(c(1,2),1)
        if(j == 1){
          substr(s,offset,offset+nchar(add.sequence[[i]])) <- add.sequence[[i]] 
        }else{
          substr(s,offset,offset+nchar(add.sequence[[i]])) <- as.character(reverseComplement(DNAString(x=add.sequence[[i]])))
        }
      }
    }
  }
  
  # add a PFM match (sampled from the pfm matrix according to frequencies)
  if("matrix" %in% is(add.pfm.match) | "data.frame" %in% is(add.pfm.match) | "list" %in% is(add.pfm.match)){
    
    if(!"list" %in% is(add.pfm.match)){  # make into temp list if single PFM
      add.pfm.match <- list(add.pfm.match)
    }
      
    for(i in c(1:length(add.pfm.match))){  # add a sequence for every entrie in the list
    
      if("matrix" %in% is(add.pfm.match[[i]])){
        add.pfm.match[[i]] <- as.data.frame(add.pfm.match[[i]])
      }
      
      # sample seq from PFM
      sampled.match <- paste0(apply(add.pfm.match[[i]], 2, function(x){
          s <- sample(c("A","C","G","T"), 1, replace=TRUE, prob=x)
        }), collapse="")
      
      offset <- sample(c(1:(l-nchar(sampled.match))), 1) # sample offset where to add
      # add
      if(!add.reverse){
        # add straight forward
        substr(s,offset,offset+nchar(sampled.match)) <- sampled.match  
      }else if(add.reverse){
        #randomly decide if to add forward or reverse complement
        j <- sample(c(1,2),1)
        if(j == 1){
          substr(s,offset,offset+nchar(sampled.match)) <- sampled.match    
        }else{
          substr(s,offset,offset+nchar(sampled.match)) <- as.character(reverseComplement(DNAString(x=sampled.match)))
        }
      }
    }
  }  
  return(s)
}

