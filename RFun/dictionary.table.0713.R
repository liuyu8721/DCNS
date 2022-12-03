dictionary.table <- function(dictionary.tree) {
  
  dictionary <- data.frame()
  paths.from.Origin <- shortest_paths(dictionary.tree, from = "Origin")$vpath
  for (i in 1:length(paths.from.Origin)) {
    x <- igraph::as_ids(paths.from.Origin[[i]])
    dictionary <- rbind(dictionary,
                        data.frame(vertex.id = i,
                                   depth = paste0("community.", 1:length(x) - 1),
                                   community = x))
  }
  dictionary <- dictionary %>%
    pivot_wider(names_from = depth, values_from = community)
  dictionary <- dictionary[with(dictionary, order(community.1, community.2, community.3, 
                                                  community.4, community.5, community.6, 
                                                  community.7, community.8, community.9,
                                                  community.10, community.11, community.12,
                                                  na.last = F)), -1]
  
  return(dictionary)
}