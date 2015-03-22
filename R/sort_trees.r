## Takes a set of unique trees and calculates the weights. Then sorts them and  adds trees to 'sorted' until the sum of weights  <= alpha
sort_trees <- function(utrees, N = 4096, alpha = .95){ 
  utrees.trees <- utrees[[1]]
  hits <- as.vector(utrees[[2]]) # "hit" information from unique_trees()
  s.hits <- sort(hits, decreasing = TRUE)
  sweights <- s.hits/sum(s.hits)
  n <- length(hits)
  #  tab <- data.frame( name = names(utrees.trees), freq = hits)
  pos <- rep(NA, n)
  for(j in 1:n){
    all.ps <- grep(s.hits[j], hits)
    u.ps <- setdiff(all.ps, na.omit(pos))
    if(length(u.ps)==1){
      pos[j] <- u.ps 
    } else{
      random.p <- sample(u.ps, 1) # breaking ties randomly
      pos[j] <- random.p
      
    }  
  }
  sorted <- list() # hate objects of varying size, but...
  class(sorted) <- class(utrees.trees)  
  i <- 1 
  sum <- 0
  while(sum <= alpha && i <= N){
    sum <- sum + as.numeric(sweights[i])
    sorted[[i]] <- utrees.trees[[pos[i]]]
    i <- i + 1   
  }
  cat("Total mass: ", sum, "\n")
  names(sorted) <- names(utrees.trees)[pos[1:(i-1)]]
  return(sorted)
}
