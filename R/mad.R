#' mad
#'
#' @param M cell-type proportions
#' @param M_est deconvolution results
#'
#'@return mean absolute difference (MAD) between true and estimated proportions
#'
#'
#'
#'@export

mad = function(M,M_est){
  return(mean(abs(M-M_est)))
}
