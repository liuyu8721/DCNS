rm(list=ls())
library(tidyr)
library(tidyverse)
library(CooccurrenceAffinity)

################################################################################

setwd("...")
source("./RFun/affinity.def0703.R")
load("./RData/calendar.RData")

################################################################################

for (w in calendar$colWeek) {

  ################################################################################
  # metadata
  ################################################################################
  load(paste0("./RData/ByWeek/w.meta.", w, ".RData"))
  print(table(w.meta$Clade))

  ################################################################################
  # Co-occurrence matrix
  ################################################################################
  load(paste0("./cooccur/cooccur.matrix.", w, ".RData"))
  muta.sel <- which(diag(as.matrix(cooccur.matrix)) / nrow(w.meta) > 0.01)
  cooccur.matrix <- cooccur.matrix[names(muta.sel), names(muta.sel)]

  ################################################################################
  # The affinity model
  # require(devtools)
  # install_github("kpmainali/CooccurrenceAffinity")
  ################################################################################
  affinity.est <- affinity.def(occur.mat = cooccur.matrix, N = nrow(w.meta), lev = 0.9)
  affinity.est.stat <- affinity.est$all %>%
    mutate(entity_1_prev_mA = entity_1_count_mA / nrow(w.meta) * 100,
           entity_2_prev_mB = entity_2_count_mB / nrow(w.meta) * 100,
           q_value = p.adjust(p_value, method = "BH")) %>%
    mutate(sig = (q_value < 0.001))

  save(affinity.est.stat, file = paste0("./AffinityModel/affinity.est", w, ".RData"))
}

