#' prop_rdirichlet
#'
#' @param seed seed set
#' @param n number of samples
#' @param k number of cell types
#' @param types cell types
#'
#'@import gtools
#'
#'@export


prop_rdirichlet = function(seed,n,k,types){
  set.seed(seed)
  M = t(rdirichlet(n, rep(1,k)))
  rownames(M) = types
  return(M)
}


