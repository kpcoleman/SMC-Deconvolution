#' scatter_plot
#'
#' @param df deconvolution results (columns = 'true','est','type','color')
#'
#'
#'
#'@export



scatter_plot = function(df){
  plot(df$true, df$est,col = df$color, pch = 16, cex = 2, ylim = c(0,1), xlab = 'True Mixture Proportion', ylab = 'Estimated mixture proportion')
  lines(seq(0,1,0.01),seq(0,1,0.01), lty = 3)
  legend(0.8,0.2, legend = c('ideal',rownames(M)[1],rownames(M)[2]), col = c('black','magenta','green'), lty = c(3,NA,NA), pch = c(NA,16,16))
  legend(0.1,0.9,legend = parse(text=sprintf('paste(r,\' = %s\')',round(cor(df$true, df$est),2))))
}
