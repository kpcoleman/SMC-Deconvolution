#' scatter_plot
#'
#' @param M cell-type proportions
#' @param M_est deconvolution results
#'
#'
#'@import graphics
#'@import stats
#'
#'@export



scatter_plot = function(M,M_est){
  M_combined = c()
  M_est_combined = c()
  type_combined = c()
  for (k in 1:nrow(M)){
    M_combined = append(M_combined,M[k,])
    M_est_combined = append(M_est_combined,M_est[k,])
    type_combined = append(type_combined,rep(rownames(M)[k], ncol(M)))
  }
  df = data.frame(M_combined)
  colnames(df) = 'true'
  df$est = M_est_combined
  df$type = type_combined
  if (nrow(M)==2){
    df$color = c(rep('magenta',ncol(M)),rep('green',ncol(M)))
    plot(df$true, df$est,col = df$color, pch = 16, cex = 2, xlim = c(0,1), ylim = c(0,1), xlab = 'True Mixture Proportion', ylab = 'Estimated mixture proportion')
    lines(seq(0,1,0.01),seq(0,1,0.01), lty = 3)
    legend(0.8,0.2, legend = c('ideal',rownames(M)[1],rownames(M)[2]), col = c('black','magenta','green'), lty = c(3,NA,NA), pch = c(NA,16,16))
    legend(0.1,0.9,legend = parse(text=sprintf('paste(r,\' = %s\')',round(cor(df$true, df$est),2))))
  }
  else{
    plot(df$true, df$est,col = as.factor(df$type), pch = 16, cex = 2, xlim = c(0,1), ylim = c(0,1), xlab = 'True Mixture Proportion', ylab = 'Estimated mixture proportion')
    lines(seq(0,1,0.01),seq(0,1,0.01), lty = 3)
    legend(0.8,0.2, legend = c('ideal',unique(type_combined)), col = c('black',as.factor(unique(type_combined))), lty = c(3,rep(NA,length(unique(type_combined)))), pch = c(NA,rep(16,unique(type_combined))))
    legend(0.1,0.9,legend = parse(text=sprintf('paste(r,\' = %s\')',round(cor(df$true, df$est),2))))
  }

}
