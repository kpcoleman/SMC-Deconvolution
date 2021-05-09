#' combine_est_props
#'
#' @param M_est a list containing matrices with deconvolution results
#' @param M a matrix of cell-type proportions
#'
#'
#'@import graphics
#'@import stats
#'
#'@export

combine_est_props = function(M_est,M){
  f1 = function(M_est1,M1){
    diff1 = sum(abs(M1-M_est1))
    diff2 = sum(abs(M1-matrix(c(M_est1[2,],M_est1[1,]), nrow = 2, byrow = T)))
    if (diff1<=diff2){
      return(M_est1)
    }
    else{
      return(matrix(c(M_est1[2,],M_est1[1,]), nrow = 2, byrow = T))
    }
  }
  M_est2 = lapply(M_est,f1,M)
  M_est2 = Reduce('+', M_est2)
  return(M_est2/length(M_est))
}
