# Takes a set of unique trees and calculates the weights. Then sorts them and  adds trees to 'sorted' until the sum of weights  <= alpha
sort_trees <- function(trees, N = 4096, alpha = .95){ 
 hits <- as.vector(trees[[2]]) # "hit" information from unique_trees()
 trees <- trees[[1]]
 weights <- hits/sum(hits) # normalised weights 
 sweights <- sort(weights, decreasing = TRUE)
 pos <- sapply(sweights, function(x) match(x, weights)) 
   #
 sorted <- list() # hate objects of varying size, but...
 class(sorted) <- class(trees)
 i <- 1 
 sum <- 0
 # 
 while(sum <= alpha && i <= N){
   sum <- sum + as.numeric(sweights[i])
   sorted[[i]] <- trees[[pos[i]]]
     i <- i + 1   
 }
 cat("Total mass: ", sum, "\n")
 return(sorted)
}
