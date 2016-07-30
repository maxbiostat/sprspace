library(ape)
library(phybase)
#################
nw.order <- function(tree){
  require(ape)
  ## please note that it depends on nw_order, and assumes it lives in /usr/bin/nw_order
  tmp <- tempfile("ordered_trees", fileext = ".nwk")
  system(input = write.tree(tree), command = paste("nw_order - >", tmp) )
  res <- read.tree(tmp)
  return(res)
}
# nw.order(rmtree(n = 50, N = 100))

tree2string <- function(tree) {tree$edge.length <- NULL; nw <- write.tree(tree); return( gsub(":[0-9]+", "", gsub("\\.", "", nw)))} 

unique_trees <- function(trees){
  urep.trees <- nw.order(trees)
  tt <- table(as.character(lapply(urep.trees, tree2string)))
  new.trees <- lapply(names(tt), function(t) read.tree(text = t))
  class(new.trees) <- "multiPhylo"
  res <- list(trees = new.trees,  hits  = as.numeric(tt))
  return(res)
}
#################
## Testing
# N <- 15
# t <- 500 
# M <- 100
# utrees <- rmtree(n = t, N = N, rooted = TRUE)
# fake.posterior <- vector(M, mode = "list")
# class(fake.posterior) <- class(utrees)
# for(i in 1:M){
#   fake.posterior[[i]] <- utrees[[sample(1:N, 1)]]  
# }
# write.tree(fake.posterior, file = "~/test_posterior_1494.nwk")
# system.time(
#   Utrees <- unique_trees(fake.posterior)
# )
# 
# ## real test
# source("../../GREPO/R/read.nexusB.r")
# real.posterior <- read.nexusB("../../DATA/dengue/Dengue4.env.trees")
# system.time(
#   uposterior <- unique_trees(real.posterior)
#   )
# max(uposterior[[2]])
