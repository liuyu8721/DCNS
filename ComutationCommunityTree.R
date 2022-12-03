rm(list=ls())
library(tidyverse)
library(igraph)
library(stringi)

setwd("...")
source("./RFun/community.network.develop.0531.R")
source("./RFun/dictionary.table.0713.R")
source("./RFun/community.similarity.R")
source("./RFun/community.refinement.0713.R")

################################################################################
# Load the current co-mutation community list
################################################################################
w <- "2020-10"
fpath <- paste0("./data/CoNet/w.", w, ".CoNet.list.csv")
community.list <- read.csv(fpath, header = T)
community.vector <- community.list$community

#-------------------------------------------------------------------------------
# Weekly co-mutation community tree
#-------------------------------------------------------------------------------
community.network <- community.network.develop(w = w, community.vector,
                                               included.cutoff =  0.88, setdiff.cutoff = 0.88)
plot(community.network,
     layout = layout_as_tree(community.network, flip.y = T),
     vertex.label = NA,
     vertex.label.cex = 0.6,
     vertext.size = 0.1,
     arrow.size = 0.5,
     arrow.width = 0.5)
