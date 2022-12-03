rm(list=ls())
library(tidyr)
library(tidyverse)

################################################################################

setwd("...")
load("./RData/calendar.RData")

################################################################################

muta.list <- vector(mode = 'list', length = nrow(calendar))
names(muta.list) <- calendar$colWeek
for (w in calendar$colWeek) {
  
  print(w)
  ################################################################################
  # read in text and deduplicate
  ################################################################################
  load(paste0("./RData/ByWeek/w.meta.", w, ".RData"))
  print(table(w.meta$Clade))
  
  load(paste0("./RData/ByWeek/w.muta.", w, ".RData"))
  length(unique(w.muta$sample))
  
  w.muta.modify <- w.muta %>%
    mutate(mutation = paste0(refvar, refpos, qvar),
           TF = 1) %>%
    select(sample, TF, mutation) 

  w.muta.modify2=as.data.frame(w.muta.modify %>% group_by(mutation) %>% count())
  w.muta.modify3=w.muta.modify2[w.muta.modify2$n>nrow(w.meta)*0.01,]
  # w.muta.modify3=w.muta.modify2
  dim(w.muta.modify3)
  w.muta.modify4=w.muta.modify[w.muta.modify$mutation %in% w.muta.modify3$mutation,]
  muta.matrix <- w.muta.modify4 %>% 
    pivot_wider(names_from = mutation, values_from = TF, values_fn = length, values_fill = 0)
  # save(muta.matrix, file = paste0("./mutation/muta.matrix.", w, ".RData"))

  L <- ncol(muta.matrix) - 1
  mutation <- colnames(muta.matrix)[-1]
  cooccur.matrix <- matrix(NA, nrow = L, ncol = L, dimnames = list(mutation, mutation))
  for (i in 1:(L-1)) {
    for (j in (i+1):L) {
      cooccur.matrix[mutation[i], mutation[i]] <- sum(muta.matrix[, mutation[i]]==1)
      
      cooccur <- muta.matrix[, mutation[i]] + muta.matrix[,mutation[j]]
      cooccur.matrix[mutation[i], mutation[j]] <- cooccur.matrix[mutation[j], mutation[i]] <- sum(cooccur==2)
    }
  }
  cooccur.matrix[mutation[L], mutation[L]] <- sum(muta.matrix[, mutation[L]]==1)
  # save(cooccur.matrix, file = paste0("./cooccur/cooccur.matrix.", w, ".RData"))
  
  
  rm(w.meta)
  rm(w.muta.modify)
  rm(w.muta)
}

