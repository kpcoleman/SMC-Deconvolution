#' jsd_mse
#'
#' @param M cell-type proportions
#' @param M_est deconvolution results
#'
#'@return matrix containing JSD between true and estimated proportions for each cell type.
#'
#'@import philentropy
#'@import MLmetrics
#'
#'@export

jsd_mse = function(M,M_est){
  jsd1 = JSD(matrix(c(M[1,],M_est[1,]), nrow = 2, byrow = T))
  jsd2 = JSD(matrix(c(M[2,],M_est[2,]), nrow = 2, byrow = T))
  mse1 = MSE(M[1,],M_est[1,])
  mse2 = MSE(M[2,],M_est[2,])
  eval = matrix(c(jsd1, jsd2, mse1, mse2), nrow = 2)
  rownames(eval) = rownames(M)
  colnames(eval) = c('JSD', 'MSE')
  return(eval)
  }
