#' mad_plot
#'
#' @param mad vector of MAD values
#' @param size vector of sample size
#'
#'@return plot of MAD between true and estimated proportions for different sample sizes
#'
#'
#'
#'@export

mad_plot = function(mad,size){
  return(barplot(mad, xlab = 'Sample size', ylab = 'Average MAD (average cell type proportions)', names.arg = size, col = 'blue'))
}
