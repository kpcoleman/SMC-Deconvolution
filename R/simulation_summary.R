#' simulation_summary
#'
#' @param t a vector of run times
#'
#'@return mean of t
#'@return standard deviation of t
#'
#'
#'@import stats
#'
#'@export

simulation_summary = function(t){
  return(c(mean(t), sd(t)))
}
