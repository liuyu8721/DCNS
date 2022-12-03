environment(CooccurrenceAffinity::affinity)

affinity.def <- function (occur.mat, N, sigdigit = NULL, lev = 0.95, ...) {

  sum(1:(ncol(occur.mat) - 1))
  output.long <- data.frame(matrix(ncol = 17, nrow = sum(1:(ncol(occur.mat) - 1))))
  colnames(output.long) <- c("entity_1", "entity_2", 
                             "entity_1_count_mA", "entity_2_count_mB", 
                             "obs_cooccur_X", "total_N", "p_value", 
                             "exp_cooccur", "alpha_mle", 
                             "conf_level", "ci_blaker", "scal", 
                             "ochiai", "jaccard", "sorensen", 
                             "simpson", "errornote")
  
  for (col_main in 1:(ncol(occur.mat) - 1)) {
    for (col_pair in (col_main + 1):ncol(occur.mat)) {
     
      X <- occur.mat[col_pair, col_main]
      mA <- occur.mat[col_main, col_main]
      mB <- occur.mat[col_pair, col_pair]

      serialnum <- (col_main - 1) * ncol(occur.mat) + col_pair - sum(1:col_main)
      if (is.null(sigdigit)) sigdigit <- 3
      
      #out.column.1 ~ out.column.6
      output.long$entity_1[serialnum] <- colnames(occur.mat)[col_main]
      output.long$entity_2[serialnum] <- colnames(occur.mat)[col_pair]
      output.long$entity_1_count_mA[serialnum] <- sum(mA)
      output.long$entity_2_count_mB[serialnum] <- sum(mB)
      output.long$obs_cooccur_X[serialnum] <- X
      output.long$total_N[serialnum] <- N
      
      #out.column.13 ~ out.column.16
      myochiai <- round(X/sqrt(as.numeric(sum(mA))*as.numeric(sum(mB))), sigdigit)
      myjaccard <- round(X/(sum(mA) + sum(mB) - X), sigdigit)
      mysorensen <- round(2 * X/(sum(mA) + sum(mB)), sigdigit)
      mysimpson <- round(X/min(sum(mA), sum(mB)), sigdigit)
      output.long$ochiai[serialnum] <- myochiai
      output.long$jaccard[serialnum] <- myjaccard
      output.long$sorensen[serialnum] <- mysorensen
      output.long$simpson[serialnum] <- mysimpson
      
      ##########################################################################
      # Maximum likelihood estimate and intervals of alpha <- ML.Alpha
      ##########################################################################
      
      ML.Alpha.def <- function (x, marg, lev = 0.95) {
        require(BiasedUrn)
        mA = marg[1]
        mB = marg[2]
        N = marg[3]
        xmin = max(mA + mB - N, 0)
        xmax = min(mA, mB)
        
        if (x < 0 | x < mA + mB - N | x > min(mA, mB)) return("Impossible x!")
        if (length(intersect(c(mA, mB), c(0, N)))) return("Degenerate co-occurrence distribution!")
        
        # scal = min(log(2 * marg[3]^2), 10)
        scal = log(2 * marg[3]^2)

        # The consistent estimator of alpha
        if (x == xmin) {
          est <- -scal 
          CI.lower.exp <- NA
          CI.upper.exp <- NA
        } else if (x == xmax) {
          est <- scal
          CI.lower.exp <- NA
          CI.upper.exp <- NA
        } else {
          # The consistent estimator for alpha
          P1 <- x/mA
          P2 <- (mB-x)/(N-mA)
          t <- (P1*(1-P2))/(P2*(1-P1))
          est <- log(t)
          
          # An approximate (1 - level) confidence interval for exp(alpha)
          H <- sqrt(1/(mA*P1*(1-P1)) + 1/((N-mA)*P2*(1-P2)))
          Z.lev <- qnorm((1+lev)/2)
          CI.lower.exp <- ((x-Z.lev/H)*(N-mA-mB+x-Z.lev/H)) / ((mA-x+Z.lev/H)*(mB-x+Z.lev/H))
          CI.upper.exp <- ((x+Z.lev/H)*(N-mA-mB+x+Z.lev/H)) / ((mA-x-Z.lev/H)*(mB-x-Z.lev/H))
        }

         list(est = est, 
             lev = lev, 
             CI.Blaker = c(CI.lower.exp, CI.upper.exp), 
             scal = scal,
             Null.Exp = as.numeric(mA) * as.numeric(mB) / N, 
             pval = AcceptAffin(x, marg, 0))
      }
      
      marg = c(mA, mB, N)
      mle <- ML.Alpha.def(x = X, marg, lev = lev)
      if (mle[1] == "Degenerate co-occurrence distribution!") {
        output.long$errornote[serialnum] <- "Degenerate co-occurrence distribution!"
        next
      }
      else {
        if (mle$pval > 1e-04) {
          output.long$p_value[serialnum] <- round(mle$pval, 4)
        }
        if (mle$pval <= 1e-04) {
          output.long$p_value[serialnum] <- formatC(mle$pval, format = "e", digits = 4)
        }
        output.long$exp_cooccur[serialnum] <- round(mle$Null.Exp, sigdigit)
        output.long$alpha_mle[serialnum] <- round(mle$est, sigdigit)
        output.long$conf_level[serialnum] <- mle$lev
        output.long$ci_blaker[serialnum] <- paste0("[", round(mle$CI.Blaker[1], sigdigit), ", ", round(mle$CI.Blaker[2], sigdigit), "]")
        output.long$scal[serialnum] <- mle$scal
      }
    }
  }
  finalout <- list(all = output.long, occur_mat = occur.mat)

  message("~~~~~~~~~~ printing head of all elements of the output list ~~~~~~~~~~")
  print(sapply(finalout, head))
  message("~~~~~~~~~~ COMPLETED: printing head of all elements of the output list ~~~~~~~~~~")
  invisible(return(finalout))
}
