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
# Load the latest co-mutation community dictionary
################################################################################
w <- "2020-09"
dictionary.old.table <- read.csv(paste0("./dictionary/dictionary.", w, ".csv"), header = T)

# Phylogenetic relationship
dictionary.old <- data.frame()
for (i in 1:nrow(dictionary.old.table)) {
  x <- stri_remove_empty(as.character(dictionary.old.table[i,]))
  if (length(x) > 1) {
    for (j in 2:length(x)) {
      dictionary.old <- rbind(dictionary.old,
                              data.frame(from = x[j-1],
                                         to = x[j],
                                         es = paste(x[j-1], ">", x[j])))
    }
  }
}
dictionary.old <- dictionary.old[!duplicated(dictionary.old$es),] %>%
  select(from, to)
dictionary.old <- graph_from_data_frame(dictionary.old, directed = T)

# Record V & E
dictionary.old.vertex <- igraph::as_ids(V(dictionary.old))
dictionary.old.edge <- igraph::as_data_frame(dictionary.old)
dictionary.old.edge2 <- with(dictionary.old.edge, paste0(from, ">", to))

################################################################################
# Load the current co-mutation community list
################################################################################
w <- "2020-10"
fpath <- paste0("./data/CoNet/w.", w, ".CoNet.list.csv")
community.list <- read.csv(fpath, header = T)
community.vector <- community.list$community

# Vertexes in dictionary
# community.corpora <- union(dictionary.old.vertex, community.list$community)

################################################################################
# Weekly co-mutation community network surveillance
################################################################################

#-------------------------------------------------------------------------------
# Pre-processing new communities relative to the old dictionary
# Condition: compression, extension, intersection
#-------------------------------------------------------------------------------
community.old <- setdiff(dictionary.old.vertex, community.vector)
community.new <- setdiff(community.vector, dictionary.old.vertex)
community.new
if (length(community.new) > 0) {
  
  # Community similarity between new weekly data and up-to-date dictionary
  community.likeness <- community.similarity(community.new, community.old)
  community.likeness <- community.likeness[with(community.likeness, order(community.new, -Jaccard)),]
  
  if (nrow(community.likeness)) {
    
    community.likeness$community.update <- NA
    for (cm in unique(community.likeness$community.new)) {
      
      i <- with(community.likeness, which(community.new==cm))
      input <- community.likeness[i,]
      community.likeness[i, "community.update"][1] <- community.refinement(input)

    }
    
    # Update community.corpora
    community.vector <- union(setdiff(community.vector, community.likeness$community.new), 
                              stri_remove_empty(str_split(na.omit(community.likeness$community.update), ";", simplify = T)))
  }
} else {
  community.likeness <- data.frame()
}
write.csv(community.likeness, file = paste0("D:/00SARS-Cov-2_RCM/ADS/Overlap/", w, ".community.overlap.csv"), row.names = F, na = "")

#-------------------------------------------------------------------------------
# Include all paternal communities relative to the old dictionary
#-------------------------------------------------------------------------------
community.vector.update <- community.vector
for (i in 1:length(community.vector)) {
  loc <- which(dictionary.old.table==community.vector[i], arr.ind=TRUE)
  if (nrow(loc)) {
    community.vector.update <- union(community.vector.update, as.character(dictionary.old.table[loc[1,1], loc[1,2]:1]))
  }
}

#-------------------------------------------------------------------------------
# Weekly co-mutation community network
#-------------------------------------------------------------------------------
community.network <- community.network.develop(w = w, community.vector = setdiff(community.vector.update, c("Origin")),
                                               included.cutoff =  0.88, setdiff.cutoff = 0.88)
community.network.vertex <- igraph::as_ids(V(community.network)) 
community.network.edge <- igraph::as_data_frame(community.network)
community.network.edge2 <- with(community.network.edge, paste0(from, ">", to))

community.network.vertex.label <- community.network.vertex
if (!setequal(intersect(dictionary.old.vertex, community.network.vertex.label), community.network.vertex.label))
  community.network.vertex.label[community.network.vertex %in% 
                                  setdiff(dictionary.old.vertex, str_split(community.network.edge2[!(community.network.edge2 %in% dictionary.old.edge2)], ">", simplify = T)[,2])] <- NA
plot(community.network,
     vertex.color = factor(community.network.vertex %in% dictionary.old.vertex, 
                           levels = c("TRUE", "FALSE")),
     edge.color = factor(community.network.edge2 %in% dictionary.old.edge2, levels = c("TRUE", "FALSE")),
     vertex.label = NA,
     vertex.label.cex = 0.6,
     vertext.size = 0.1,
     arrow.size = 0.5,
     arrow.width = 0.5)

plot(community.network,
     layout = layout_as_tree(community.network, flip.y = T),
     vertex.label = community.network.vertex.label,
     vertex.color = factor(community.network.vertex %in% dictionary.old.vertex,
                           levels = c("TRUE", "FALSE")),
     edge.color = factor(community.network.edge2 %in% dictionary.old.edge2, levels = c("TRUE", "FALSE")),
     vertex.label.cex = 0.8)

################################################################################
# Update the dynamic dictionary
################################################################################

#-------------------------------------------------------------------------------
# Phylogenetic relationship comparison between dictionary tree and current network
#-------------------------------------------------------------------------------
# community.network.table <- dictionary.table(community.network)

# Search contradictory phylogenetic relationship
community.old <- intersect(community.network.vertex, dictionary.old.vertex)
community.paths <- shortest_paths(community.network, from = "Origin", to = community.old)$vpath
dictionary.paths <- shortest_paths(dictionary.old, from = "Origin", to = community.old)$vpath

community.paths.masked <- data.frame()
for (i in 1:length(community.paths)) {
  x <- paste0(as_ids(community.paths[[i]]), collapse = ">")
  y <- paste0(as_ids(dictionary.paths[[i]]), collapse = ">")
  if (x!=y) {
    community.paths.masked <- rbind(community.paths.masked,
                                    data.frame(new = x,
                                               old = y,
                                               update = y))
  }
}

# Update the community network
if (nrow(community.paths.masked)) {
  for (j in 1:nrow(community.paths.masked)) {
    #delete edges
    y <- str_split(community.paths.masked$new[j], pattern = ">")[[1]]
    for (i in 1:(length(y) - 1)) {
      community.network.edge <- community.network.edge %>%
        filter(!(from==y[i] & to==y[i+1]))
    }
    #add edges
    x <- str_split(community.paths.masked$update[j], pattern = ">")[[1]]
    for (i in 1:(length(x) - 1)) {
      community.network.edge <- rbind(community.network.edge,
                                      data.frame(from = x[i],
                                                 to = x[i+1]))
    }
  }
}
community.network.edge <- community.network.edge %>%
  mutate(edge = paste0(from, ">", to)) %>%
  filter(!duplicated(edge)) %>%
  select(from, to)

community.network <- igraph::graph_from_data_frame(community.network.edge)
plot(community.network,
     layout = layout_as_tree(community.network, flip.y = T),
     vertex.label = NA)

#-------------------------------------------------------------------------------
# Update the dictionary
#-------------------------------------------------------------------------------
dictionary.tree <- igraph::union(community.network, dictionary.old)
plot(dictionary.tree, 
     layout = layout_as_tree(dictionary.tree, flip.y = T),
     vertex.label = NA)

# Tree update
dictionary.vertex <- igraph::as_ids(V(dictionary.tree))
dictionary.edge <- igraph::as_data_frame(dictionary.tree)
dictionary.edge.update <- (with(dictionary.edge, paste0(from, ">", to)) %in% dictionary.old.edge2)
dictionary.vertex.update <- (dictionary.vertex %in% dictionary.old.vertex)

dictionary.vertex.label <- dictionary.vertex
dictionary.vertex.label[dictionary.vertex %in% dictionary.old.vertex] <- NA

plot(dictionary.tree, 
     vertex.color = factor(dictionary.vertex.update, levels = c("TRUE", "FALSE")),
     edge.color = factor(dictionary.edge.update, levels = c("TRUE", "FALSE")),
     layout = layout_as_tree(dictionary.tree, flip.y = T),
     vertex.label = dictionary.vertex.label,
     vertex.label.cex = 0.8,
     main = w)

save(dictionary.tree, file = paste0("./data/dictionary/dictionary.", w, ".RData"))

dictionary.new.table <- as.matrix(dictionary.table(dictionary.tree))
dictionary.new.table <- dictionary.new.table[with(as.data.frame(dictionary.new.table), order(community.1, community.2, community.3, community.4, community.5, na.last = F)),]
write.csv(dictionary.new.table, file = paste0("./data/dictionary/dictionary.", w, ".csv"), 
          quote = T, row.names = F, na = "")

