#' smc_deconvolution
#'
#' @param Y bulk GE data
#' @param K number of cell types
#' @param N number of samples
#' @param L number of iterations
#' @param parallel whether or not to use parallel computing
#' @param num_cores number of cores
#'
#'@return a data frame with deconvolution results
#'@return runtime of method
#'
#'
#'@importFrom stats rgamma rnorm runif
#'@importFrom pracma randsample
#'@importFrom rlang duplicate
#'@import parallel
#'@import foreach
#'@import doParallel
#'
#'@export


smc_deconvolution = function(Y,K,N=40,L=1, parallel= F, num_cores=1){
  if (K<2){
    stop('Number of cell types must be >= 2')
  }
  if (num_cores<1 || num_cores > detectCores()){
    stop(paste('Number of cores must be between 1 and'),detectCores())
  }
  if (N<1){
    stop('Number of samples must be >= 1')
  }
  if (L<1){
    stop('Number of iterations must be >= 1')
  }
  deconvolution = function(Y,K,N){
    J=ncol(Y)
    I = nrow(Y)
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
  if (parallel==T){
    start.time = Sys.time()
    cl = makeCluster(num_cores,type="SOCK")
    registerDoParallel(cl)
    props = foreach(l=1:L) %dopar% {
      deconvolution(Y,K,N)
    }
    stopCluster(cl)
    end.time = Sys.time()
    runtime = end.time - start.time
  }
  else {
    props = vector(mode = "list", length = L)
    runtime = c(rep(NA,L))
    for (l in 1:L){
      print(paste('Itertation',l))
      start.time = Sys.time()
      props[[l]] = deconvolution(Y,K,N)
      end.time = Sys.time()
      runtime[l] = end.time - start.time
    }
  }
  #props_est = Reduce('+',props)/L
  return(list(props,runtime))
}







