#' prop_rdirichlet
#'
#' @param seed seed set
#' @param n number of samples
#' @param K number of cell types
#' @param types cell types
#'
#'@import gtools
#'
#'@export


prop_rdirichlet = function(seed,n,K,types){
  set.seed(seed)
  M = t(rdirichlet(n, rep(1,K)))
  rownames(M) = types
  return(M)
}


