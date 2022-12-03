community.compression <- function(c.new, c.dict) {
  return(c.dict)
}
#----------- Example 1
# community.compression("A17858G|C18060T",
#                       "C17747T|A17858G|C18060T")
# > "C17747T|A17858G|C18060T"

community.extension <- function(c.new, c.dict) {
  c.new.muta <- c(str_split(c.new, "\\|", simplify = T))
  c.dict.muta <- c(str_split(c.dict, "\\|", simplify = T))
  c.diff.muta <- setdiff(c.new.muta, c.dict.muta)
  if (length(c.diff.muta)>1) return(paste(c.dict, paste(c.diff.muta, collapse = "|"), sep = ";"))
  else return(c.dict)
}
#----------- Example 1
# community.extension("C6312A|C13730T|C23929T",
#                     "C6312A|C23929T")
# > "C6312A|C23929T"
#----------- Example 2
# community.extension("C6312A|C13730T|C23929T|C28311T",
#                     "C6312A|C23929T")
# > "C6312A|C23929T;C13730T|C28311T"

community.refinement <- function(input) {
  
  # Simple community compression
  if (nrow(input)==1 & input$relate[1]=="A<B") {
    return(input$community.old)
  }
  
  # Simple community extension
  if (nrow(input)==1 & input$relate[1]=="A>B") {
    return(community.extension(input$community.new, input$community.old))
  }
  
  # Simple community overlap
  if (nrow(input)==1 & input$relate[1]=="AnB") {
    if (input$Jaccard>=0.5) return(input$community.old) else return(input$community.new)
  }
  
  if (nrow(input) > 1) {
    
    c.update <- c()
    
    # priority for extension
    input.extent <- subset(input, relate=="A>B")
    if (nrow(input.extent)) {
      for (i in 1:nrow(input.extent)) {
        c.new <- input.extent$community.new[i]
        c.dict <- input.extent$community.old[i]
        c.solution <- community.extension(c.new, c.dict)
        c.update <- c(c.update, word(c.solution, 1, 1, fixed(";")))
        if (i < nrow(input.extent)) input.extent$community.new[i+1] <- word(c.solution, 2, 2, fixed(";"))
      }
      
      #update community.new
      c.rest <- word(c.solution, 2, 2, fixed(";"))
      if (!is.na(c.rest)) {
        if (nrow(subset(input, relate!="A>B"))) input$community.new[input$relate!="A>B"] <- c.rest else c.update <- c(c.update, c.rest)
      } else input <- subset(input, relate=="A>B")
    }
    
    # other condiction
    input.other <- subset(input, relate!="A>B")
    if (nrow(input.other)) {
      community.jaccard <- community.similarity(unique(input.other$community.new), input.other$community.old)
      if (nrow(community.jaccard)==0) c.update <- c(c.update, input.other$community.new)
      else if (nrow(community.jaccard) & length(intersect(community.jaccard$relate, "A<B"))) c.update <- c(c.update, input.other$community.old[community.jaccard$relate=="A<B"])
      else if (nrow(community.jaccard) & max(community.jaccard$Jaccard) < 0.5) c.update <- c(c.update, unique(input.other$community.new))
      else if (nrow(community.jaccard) & max(community.jaccard$Jaccard) >= 0.5) c.update <- c(c.update, community.jaccard$community.old[community.jaccard$Jaccard==max(community.jaccard$Jaccard)])
    }

    return(paste(c.update, collapse = ";"))
  }
}

# cm="T733C|C2749T|C3828T|A5648C|A6319G|A6613G|C12778T|C13860T|G17259T|C21614T|C21621A|C21638T|G21974T|G22132T|A22812C|C23525T|C24642T|G25088T|T26149C|G28167A|C28512G|A28877T|G28878C"
# i <- with(community.likeness, which(community.new==cm))
# input <- community.likeness[i,]
# community.refinement(input)

