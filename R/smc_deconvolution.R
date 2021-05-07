#' smc_deconvolution
#'
#' @param X cell-type specific GE data
#' @param Y bulk GE data
#' @param N number of samples
#' @param L number of iterations
#'
#'@return a data frame with deconvolution results
#'
#'
#'@importFrom stats rgamma rnorm runif
#'@importFrom pracma randsample
#'@importFrom rlang duplicate
#'
#'@export


smc_deconvolution = function(X,Y,M,N=40,L=1, num_cores=1){
  deconvolution = function(X,Y,N){
    K = ncol(X)
    J=ncol(Y)
    I = nrow(X)
    #starting values
    mu_kj = 0
    nu_kj = 1/0.01

    #mu_ik = runif(1,1,100)
    mu_ik = matrix(runif(K*I,0,100), nrow = I)
    nu_ik = 1/1000

    M1 = vector(mode = "list", length = N)
    X1 = vector(mode = "list", length = N)

    alpha = 2
    beta = 10000

    for(n in 1:N){
      X1[[n]] = matrix(rep(NA,K*I), nrow = I)
      M1[[n]] = matrix(rnorm(K*J,mu_kj,sqrt(1/nu_kj)), nrow = K)
      for (i in 1:I){
        for (k in 1:K){
          X1[[n]][i,k] = rnorm(1,mu_ik[i,k],sqrt(1/nu_ik))
        }
      }
    }

    lambda = rgamma(N,alpha,beta)
    w = rep(1/N,N)

    epsilon = seq(0,1,10e-4)
    T1 = length(epsilon)


    pdf_y = rep(NA,N)
    log_pdf_y = rep(NA,N)

    for (t in 2:T1){
      print(t)
      for (n in 1:N){
        log_pdf_y[n]= I*J*(epsilon[t]-epsilon[t-1])/2*log(lambda[n]) -lambda[n]*(epsilon[t]-epsilon[t-1])/2*sum((Y-(X1[[n]]%*%M1[[n]]))^2)
      }
      #XM = Map('%*%',X1,M1)
      #sum_Y_XM_2 = unlist(lapply(lapply(lapply(lapply(XM,'*',-1),'+',Y),'^',2),sum))
      #log_pdf_y = I*J*(epsilon[t]-epsilon[t-1])/2*log(lambda)-lambda*(epsilon[t]-epsilon[t-1])/2*sum_Y_XM_2
      pdf_y = exp(log_pdf_y-max(log_pdf_y))
      w_old = w
      w = w_old*pdf_y
      w = w/sum(w)
      ess = 1/sum(w^2)
      if (ess<N/10){
        rr = pracma::randsample(N,N,w = w, replacement = T)
        lambda = lambda[rr]
        X1 = X1[rr]
        M1 = M1[rr]
        w = rep(1/N,N)
      }
      for (n in 1:N){
        X1_n = X1[[n]]
        M1_n = M1[[n]]
        lambda_n = lambda[n]
        alpha_n = alpha+(I*J*epsilon[t])/2
        beta_n = beta+epsilon[t]*sum((Y-(X1_n%*%M1_n))^2)/2
        lambda_n =  rgamma(1,alpha_n, beta_n)
        epsilon_n = epsilon[t]*lambda_n
        lambda[n] = lambda_n
        for (k in 1:K){
          for (j in 1:J){
            X1_k = X1_n[,k]
            Y_j = Y[,j]
            U_kj = nu_kj+epsilon_n*sum(X1_k^2)
            #X2_k = duplicate(X1_n)
            X2_k = X1_n
            X2_k[,k] = 0
            M1_j = M1_n[,j]
            B_ijk = X2_k%*%M1_j
            V_kj = mu_kj*nu_kj + epsilon_n*sum(Y_j*X1_k)-epsilon_n*sum(X1_k*B_ijk)
            M1_n[k,j] = rnorm(1,V_kj/U_kj,sqrt(1/U_kj))
            if (M1_n[k,j] < 0){
              M1_n[k,j]=0
            }
          }
        }
        M1[[n]] = M1_n
        for (i in 1:I){
          for (k in 1:K){
            M1_k = M1[[n]][k,]
            Y_i = Y[i,]
            #X1_i = duplicate(X1[[n]][i,])
            X1_i = X1[[n]][i,]
            X1_i[k] = 0
            Y1_i = X1_i%*%M1[[n]]
            A_ik = nu_ik + epsilon_n*sum(M1_k^2)
            B_ik = nu_ik*mu_ik[i,k] + epsilon_n*sum(M1_k*Y_i) - epsilon_n*sum(M1_k*Y1_i)
            X1_n[i,k] = rnorm(1, B_ik/A_ik, sqrt(1/A_ik))
          }
        }
        X1[[n]] = X1_n
      }
    }
    X_est = Reduce('+', Map('*',X1,w))
    M_est1 = Reduce('+', Map('*',M1,w))
    #M_est[l] = prop.table(M_est1,2)
    M_est = prop.table(M_est1,2)
    return(M_est)
  }
  cl = makeCluster(num_cores,type="SOCK")
  registerDoParallel(cl)
  props = foreach(l=1:L) %dopar% {
    print(l)
    deconvolution(X,Y,N)
  }
  stopCluster(cl)
  props_est = Reduce('+',props)/L
  df = data.frame(c(M[1,],M[2,]))
  colnames(df) = 'true'
  df$est = c(props_est[1,],props_est[2,])
  df$type = c(rep(rownames(M)[1],ncol(M)),rep(rownames(M)[2],ncol(M)))
  df$color = c(rep('magenta',ncol(M)),rep('green',ncol(M)))
  #plot(df$true, df$est,col = df$color, pch = 16, cex = 2, ylim = c(0,1), xlab = 'True Mixture Proportion', ylab = 'Estimated mixture proportion')
  #lines(seq(0,1,0.01),seq(0,1,0.01), lty = 3)
  #legend(0.8,0.2, legend = c('ideal',rownames(M)[1],rownames(M)[2]), col = c('black','magenta','green'), lty = c(3,NA,NA), pch = c(NA,16,16))
  #legend(0.1,0.9,legend = parse(text=sprintf('paste(r,\' = %s\')',round(cor(df$true, df$est),2))))
  return(df)
}







