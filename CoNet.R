rm(list=ls())
library(tidyr)
library(tidyverse)
library(CooccurrenceAffinity)

################################################################################

setwd("...")
load("./RData/calendar.RData")

library(RColorBrewer)
color.list <- c(brewer.pal(12, "Paired"), brewer.pal(12, "Set3"))

################################################################################
for (w in calendar$colWeek) {
  load(paste0("./AffinityModel/affinity.est", w, ".RData"))
  
  est.cooccur <- affinity.est.stat %>%
    mutate(alpha.MLE = alpha_mle) %>%
    filter(q_value < 0.001 & entity_1_prev_mA > 1 & entity_2_prev_mB > 1) %>%
    # filter(obs_cooccur_X > total_N * 0.01) %>%
    filter(alpha_mle > 0) %>%
    select(entity_1, entity_2, 
           entity_1_count_mA, entity_2_count_mB, 
           entity_1_prev_mA, entity_2_prev_mB,
           obs_cooccur_X, alpha.MLE, q_value, ochiai, jaccard, sorensen, simpson)
  muta.sel <- union(est.cooccur$entity_1, est.cooccur$entity_2)
  
  ################################################################################
  
  library(igraph)
  g.df <- est.cooccur %>%
    mutate(from = entity_1, 
           to = entity_2)
  g <- graph.data.frame(g.df, directed = F) 
  plot(g, vertex.label = NA)
  
  g.v <- as_ids(V(g))
  g.e <- as_ids(E(g))
  #-------------------------------------------------------------------------------
  # plot(ecdf(est.cooccur$ochiai), col = "blue", xlim = c(0, 1), ylim = c(0, 1), main = NULL)
  # par(new = T)
  # plot(ecdf(est.cooccur$jaccard), col = "red", xlim = c(0, 1), ylim = c(0, 1), main = NULL)
  #-------------------------------------------------------------------------------
  
  ################################################################################
  # Community detection
  ################################################################################
  g.cd <- g.df %>%
    filter(ochiai > 0.9)
  g.cd <- graph.data.frame(g.cd, directed = F) 
  # plot(g.cd, vertex.label.cex = 0.8)
  
  kc <- cluster_edge_betweenness(g.cd)
  length(kc)
  sizes(kc)
  plot(kc, g.cd)
  
  suppressMessages(library(WGCNA))
  kc.colors <- labels2colors(kc$membership)
  # kc.colors <- labels2colors(kc$membership, colorSeq = color.list)
  names(kc.colors) <- kc$names
  
  g.v.color <- rep(adjustcolor("grey60", 0.5), length(g.v)) 
  names(g.v.color) <- g.v
  g.v.color[kc$names] <- kc.colors[kc$names]
  g.v.frame.color <- ifelse(g.v.color==adjustcolor("grey60", 0.5), "grey80", "red")
  
  g.e.color <- ifelse(g.df$ochiai > 0.9, "red", adjustcolor("grey88", 0.3))
  names(g.e.color) <- g.e
  for (e in g.e) {
    if (g.e.color[e]=="red") g.e.color[e] <- g.v.color[word(e, 1, 1, fixed("|"))]
  }
  
  g.v.size <- rep(1, length(g.v))
  names(g.v.size) <- g.v
  g.v.size[kc$names] <- 3
  
  setwd("D:/00SARS-Cov-2_RCM/ADS/CoNet")
  pdf(paste0("w.control.", w, ".pdf"))
  plot(kc, g.cd, col = g.v.color[kc$names],
       vertex.label = NA)
  dev.off()
  pdf(paste0("w.CoNet.", w, ".pdf"), width = 4, height = 4)
  plot(g, 
       # layout = layout_nicely(g),
       vertex.label = NA,
       vertex.color = g.v.color,
       vertex.frame.color = g.v.frame.color,
       vertex.size = g.v.size,
       edge.width = 1,
       edge.color = g.e.color)
  dev.off()
  
  #-------------------------------------------------------------------------------
  # Community-based prevalence
  #-------------------------------------------------------------------------------
  load(paste0("./RData/ByWeek/w.meta.", w, ".RData"))
  N <- nrow(w.meta)
  load(paste0("./RData/ByWeek/w.muta.", w, ".RData"))
  
  w.muta <- w.muta %>%
    mutate(mutation = paste0(refvar, refpos, qvar)) %>%
    select(sample, mutation) 
  rownames(w.muta) <- NULL
  
  muta.community.list <- data.frame()
  for (i in 1:length(kc)) {
    muta.community <- kc[[i]]
    muta.community <- muta.community[order(as.numeric(gsub('[^0-9.]', '', muta.community)))]
    # if (setequal(muta.community, c("A17858G", "C18060T"))) muta.community <- c("C17747T", "A17858G", "C18060T")
    sample.community <- c()
    for (m in muta.community) {
      s <- subset(w.muta, mutation==m) %>% pull(sample)
      if (m==muta.community[1]) sample.community <- s else sample.community <- intersect(sample.community, s)
    }
    
    meta.community <- subset(w.meta, Accession.ID %in% sample.community)
    clade.community <- data.frame(table(meta.community$Clade))
    clade.community <- clade.community[order(clade.community$Freq,decreasing = T),]
    pangolin.community <- data.frame(table(meta.community$Lineage))
    pangolin.community <- pangolin.community[order(pangolin.community$Freq, decreasing = T),]
    muta.community.list <- rbind(muta.community.list,
                                 data.frame(collection.week = sub("-","_", w),
                                            n.sample = N,
                                            modularity = modularity(kc),
                                            community = paste(muta.community, collapse = "|"),
                                            community.count = length(sample.community),
                                            community.freq = length(sample.community) / N * 100,
                                            community.clade = paste(clade.community$Var1, collapse = "|"),
                                            community.clade.count = paste(clade.community$Freq, collapse = "|"),
                                            community.pangolin = paste(pangolin.community$Var1, collapse = "|"),
                                            community.pangolin.count = paste(pangolin.community$Freq, collapse = "|")))
  }
  muta.community.list <- muta.community.list[order(muta.community.list$community.count, decreasing = T),] %>%
    filter(grepl("\\|", community))
  write.csv(muta.community.list, file = paste0("w.", w, ".CoNet.list.csv"), quote = T, row.names = F, fileEncoding = "UTF-8")
  
}


