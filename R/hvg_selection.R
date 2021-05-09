#' hvg_selection
#'
#' @param X cell-type specific GE data
#' @param p include top (100*p)% highly variable genes
#'
#'@return GE matrix with highly variable genes
#'
#'
#'@import stats
#'
#'@export


hvg_selection = function(X,p=0.05){
  compute_cv = function(x) sd(x) / mean(x)
  cv = apply(X, 1, compute_cv)
  X1 = X[rank(cv) / length(cv) > 1 - p, ]
  return(X1)
}


