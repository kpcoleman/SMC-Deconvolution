#' create_mixture_samples
#'
#' @param X cell-type specific GE data
#' @param M cell-type proportions of samples
#'
#'@return GE matrix of samples with with cell-type proportions given by M
#'
#'
#'
#'@export

create_mixture_samples = function(X,M){
  return(X%*%M)
}

