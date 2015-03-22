library(ape)
library(phybase)
#################
nw.order <- function(tree){
## please note that it depends on nw_order, and assumes it lives in /usr/bin/nw_order
  system(input = write.tree(tree), command = "nw_order -", intern = TRUE)
}
# nw.order(rmtree(n = 50, N = 100))
obtain.urep <- function(trees){
  raw.trees <- nw.order(trees)
  n <- length(raw.trees)
  utrees <- vector(n, mode = "list")
  for (i in 1:n) {utrees[[i]] <- read.tree(text = raw.trees[i])}
  class(utrees) <- "multiPhylo"
  return(utrees)
}

tree2string <- function(tree) {nw <- write.tree(tree); return( gsub(":[0-9]+", "", gsub("\\.", "", nw)))} 
unique_trees <- function(trees){
  trees <- obtain.urep(trees)
  tt <- table(as.character(lapply(trees, tree2string)))
  n <- length(tt)
  new.trees <- vector(n, mode = "list")
  for(i in 1:n){
    new.trees[[i]] <- read.tree(text = names(tt)[i])
  }
  class(new.trees) <- "multiPhylo"
  res <- list(trees = new.trees,  hits  = as.numeric(tt))
  res
}
#################
## Testing
# utrees <- rmtree(n = 17, N = 1500, rooted = TRUE)
# M <- 10000
# fake.posterior <- vector(M, mode = "list")
# class(fake.posterior) <- class(utrees)
# for(i in 1:M){
#   fake.posterior[[i]] <- utrees[[sample(1:1500, 1)]]  
# }
# write.tree(fake.posterior, file = "~/test_posterior_1494.nwk")
# system.time(
#   Utrees <- unique_trees(fake.posterior)
#   )
# 
# ## real test
# source("../../GREPO/R/read.nexusB.r")
# real.posterior <- read.nexusB("../../DATA/dengue/Dengue4.env.trees")
# system.time(
#   uposterior <- unique_trees(real.posterior)
#   )
# max(uposterior[[2]])
