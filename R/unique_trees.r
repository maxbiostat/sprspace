# Slight modification of the orignal unique.multiPhylo() to keep track of how many times it has seen a particular topology
unique_trees <- function (x, incomparables = FALSE, use.edge.length = FALSE, 
          use.tip.label = TRUE, ...) 
{
  n <- length(x)
  keep <- 1L
  old.index <- seq_len(n)
  hits <- rep(1, n)
  for (i in 2:n) {
    cat("Doing tree", i, "\n")
    already.seen <- FALSE
    for (j in keep) {
      if (all.equal(x[[j]], x[[i]], use.edge.length = use.edge.length, 
                    use.tip.label = use.tip.label)) {
        already.seen <- TRUE
        old.index[i] <- j
        hits[j] <- hits[j] + 1
        break
      }
    }
    if (!already.seen) 
      keep <- c(keep, i)
  }
  hits <- hits[keep]
  res <- x[keep]
  res <- list(res, hits)
  attr(res, "old.index") <- old.index
  res
}
