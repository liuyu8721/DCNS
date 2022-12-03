community.network.develop <- function(w, community.vector, included.cutoff =  0.88, setdiff.cutoff = 0.88) {
  
  load(paste0("D:/00SARS-Cov-2_RCM/ADS/mutation/muta.matrix.", w, ".RData"))
  muta.matrix <- data.frame(muta.matrix)
  
  community.sample <- list()
  community.nsample <- c()
  for (i in 1:length(community.vector)) {
    clique <- community.vector[i]
    # community <- strsplit(clique, "\\|")[[1]]
    # From .12345ATCG to X.12345ATCG
    community <- as.character(lapply(strsplit(clique, "\\|")[[1]], function(x){gsub("^\\.","X.",x)}))
    if (setequal(intersect(community, colnames(muta.matrix)), community)) {
      input <- muta.matrix
      for (m in community) {
        input <- input %>% 
          filter(input[,m]==1)
      }
      community.sample[[i]] <- input %>% pull(sample)
    } else {
      community.sample[i] <- list(NULL) 
    }
    names(community.sample)[i] <- clique
    community.nsample <- c(community.nsample, length(community.sample[[i]]))
  }
  community.sample <- community.sample[names(community.sample)[order(community.nsample, decreasing = T)]]
  names(community.sample)
  length(community.sample[[1]])
  
  community.relationship <- data.frame()
  community.vector <- names(community.sample)
  for (i in 1:(length(community.vector)-1)) {
    for (j in (i+1):length(community.vector)) {
      A.count = length(community.sample[[i]])
      B.count = length(community.sample[[j]])
      AB.count = length(intersect(community.sample[[i]], community.sample[[j]]))
      jaccard <- AB.count / length(union(community.sample[[i]], community.sample[[j]]))
      community.relationship <- rbind(community.relationship,
                                          data.frame(A = community.vector[i], 
                                                     B = community.vector[j],
                                                     A.count = A.count,
                                                     B.count = B.count,
                                                     AB.count = AB.count,
                                                     AB.rate = AB.count / B.count,
                                                     Jaccard = jaccard))
    }
  }
  
  # Community relationship
  community.relationship <- community.relationship %>%
    filter(AB.rate > included.cutoff & Jaccard < setdiff.cutoff)
  
  library(igraph)
  dictionary.tree <- graph_from_data_frame(community.relationship, directed = T)
  # plot(dictionary.tree, layout = layout_as_tree(dictionary.tree, flip.y = T))
  
  # dictrion tree monification
  dictionary.edge <- igraph::as_data_frame(dictionary.tree)
  dictionary.edge$tid <- NA
  dictionary.edge$hid <- NA
  for (i in 1:nrow(dictionary.edge)) {
    dictionary.edge$tid[i] <- which(community.vector==dictionary.edge$from[i])
    dictionary.edge$hid[i] <- which(community.vector==dictionary.edge$to[i])
  }
  dictionary.edge <- dictionary.edge[with(dictionary.edge, order(hid, -tid)),]
  dictionary.edge <- dictionary.edge[!duplicated(dictionary.edge$hid),]
  
  reachability.matrix <- matrix(data = NA, 
                                nrow = length(community.vector),
                                ncol = length(community.vector),
                                dimnames = list(community.vector, community.vector))
  reachability.matrix[upper.tri(reachability.matrix)] <- 0
  for (i in 1:nrow(dictionary.edge)) {
    r = dictionary.edge$from[i]
    c = dictionary.edge$to[i]
    reachability.matrix[r,c] <- 1
  }
  
  dictionary.tree.update <- graph_from_adjacency_matrix(reachability.matrix, mode = "directed", diag = F)
  
  # Add Origin edge / Level 1 nodes detection
  dist.table <- data.frame()
  for (i in 1:length(community.vector)) {
    v <- community.vector[i]
    d <- distances(dictionary.tree.update, 
                   v = i,
                   mode = "in")
    # class(d)
    dist.table <- rbind(dist.table,
                        data.frame(v = v,
                                   d = min(d[v, -which(colnames(d)==v)])))
  }
  
  # dictrion tree monification
  dictionary.edge <- rbind(data.frame(from = "Origin",
                                      to = subset(dist.table, d==Inf)$v),
                           igraph::as_data_frame(dictionary.tree.update))
  dictionary.tree.update2 <- graph_from_data_frame(dictionary.edge, directed = T)
  # plot(dictionary.tree.update2, 
  #      layout = layout_as_tree(dictionary.tree.update2, flip.y = T),
  #      vertex.label.cex = 0.7)
  
  return(dictionary.tree.update2)
}