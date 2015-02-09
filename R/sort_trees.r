# Takes a set of trees and calculates the weights. Then sorts them and  adds trees to 'sorted' until the sum of weights  <= alpha
# TODO: find if we should also put a limit on the number of trees
sort_trees <- function(trees, N = 4096, alpha = .95){ 
 posts <- abs (unlist(lapply(trees, function(tree) tree$posterior)) )
 weights <- posts/sum(posts) # normalised weights 
 sweights <- sort(weights)
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
