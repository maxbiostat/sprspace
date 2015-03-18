## Here is where we should dump all auxiliary (aka less important) functions, not to clutter the MWE script

## Compute SPR distance between the tree and the one before it in the sample
distt2 <- function(pos, x) rspr(x[[pos-1]], x[[pos]])# TODO make parallel version

## Take a "radius" matrix from rspr.matrix() and turn it into a "proper" binary incidence matrix.
binarise.matrix <- function(m){
  # takes a matrix of weights and turns into an incidence matrix
  mout <- m
  mout[.Internal(which(m > 0))] <- 1
  mout[.Internal(which(m == -1))] <- 0
  return(mout)
}