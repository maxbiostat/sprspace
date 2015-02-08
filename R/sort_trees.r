# Takes a set of trees and calculates the weights. Then sorts them and  adds trees to 'sorted' until the sum of weights  <= alpha
# TODO: find if we should also put a limit on the number of trees
sort_trees <- function(trees, N = 4096, alpha = .95){ 
 posts <- unlist(lapply(trees, function(tree) tree$posterior))
 weights <- posts/sum(posts) # normalised weights 
 sweights <- sort(weights)
 pos <- sapply(sweights, function(x) match(x, weights)) 
   #
 sorted <- list() # hate objects of varying size, but...
 class(sorted) <- class(trees)
 i <- 1 
 sum <- 0
 # && i <= N
 while(sum <= alpha){
   sum <- sum + as.numeric(sweights[i])
   sorted[i] <- trees[pos[i]]
     i <- i + 1   
 }
#  cat("Total mass: ", sum, "\n")
#   if(length(sorted)>N){sorted <- sorted[1:N]} # if we put the limit of at most N trees 
 return(sorted)
}
