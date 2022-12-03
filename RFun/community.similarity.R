community.similarity <- function(community.new, community.old, unequal = FALSE) {
  
  similarity.matrix <- matrix(data = NA, nrow = length(community.new), ncol = length(community.old), 
                              dimnames = list(community.new, community.old))
  community.relationship <- data.frame()
  for (ci in community.new) {
    for (cj in community.old) {
      mi <- str_split(ci, "\\|", simplify = T)
      mj <- str_split(cj, "\\|", simplify = T)
      
           if (setequal(mi, mj)) similarity.matrix[ci, cj] <- "A=B"                          #EQUAL
      else if (setequal(union(mi, mj), mj)) similarity.matrix[ci, cj] <- "A<B"               #SUBSET
      else if (setequal(intersect(mi, mj), mj)) similarity.matrix[ci, cj] <- "A>B"           #SUBSET
      else if (length(intersect(mi, mj))) similarity.matrix[ci, cj] <- "AnB"                 #INTERSECT

      if (!is.na(similarity.matrix[ci, cj])) community.relationship <- rbind(community.relationship,
                                                                             data.frame(community.new = ci,
                                                                                        community.old = cj,
                                                                                        relate = similarity.matrix[ci, cj],
                                                                                        Jaccard = length(intersect(mi, mj)) / length(union(mi, mj))))
    }
  }
  
  return(community.relationship)
}