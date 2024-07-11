# For this function, I assume that there is no global effect
sim_matList = function(n, rho, num_F, k_vec, num_G, F_0=TRUE){
  #browser()
  Fk = list()
  if(F_0){
    Fk[[1]] = diag(nrow=n,ncol=n)
  }
  for(k in (1:num_F)){
    X <- t(rmultinom(n, 1, rep(1/k_vec[k], k_vec[k]))) 
    Fk[[length(Fk) + 1]] = X%*%t(X)
  }
  
  Gl = list()
  Ml = list()
  Al = list()
  for(l in (1:num_G)){
    vois <- matrix(rbinom(n^2, 1, log(n)/n), n, n)
    vois <- (vois + t(vois)); diag(vois) <- 0; vois[which(vois > 1)] <- 1
    M = diag(rowSums(vois))
    A = vois
    A = (diag(1/rowSums(A))%*%A)
    Gl[[l]] = M - M%*%A*rho
    Ml[[l]] = M
    Al[[l]] = A
  }
  
  matList = list(Fk=Fk, Gl=Gl, Ml=Ml, Al=Al)
  return(matList)
}

# setting a seed for the bigger simulations
sim_cov_with_seed =function(p, Sigma,seed=NULL){
  
  n = dim(Sigma)[1]
  
  if(!is.null(seed)){
    set.seed(seed)
  }
  
  Y <- mvtnorm::rmvnorm(p, sigma=Sigma)
  corY = cor_from_standard_errors(Y)
  # corY = matrix(0,ncol=n,nrow=n)
  # for(t in 1:p){
  #   corY = corY + t(t(Y[t,]))%*%t(Y[t,])
  # }
  # corY = corY/(p-1)
  # diag(corY)=1
  # library(Matrix)
  # # 
  # # #Find nearest positive definite matrix
  # corY = nearPD(corY,corr=TRUE,conv.tol=100000)$mat
  
  #corY = cor(Y)
  #browser()
  return(list(Y=Y,corY=corY))
}

sim_cov =function(p, Sigma,seed=NULL){
  
  n = dim(Sigma)[1]
  Y <- mvtnorm::rmvnorm(p, sigma=Sigma)
  
  #browser()
  
  # Y is already standardized
  corY = cor_from_standard_errors(Y)
  
  # corY = matrix(0,ncol=n,nrow=n)
  # for(t in 1:p){
  #   corY = corY + t(t(Y[t,]))%*%t(Y[t,])
  # }
  # corY = corY/(p-1)
  # diag(corY)=1
  # library(Matrix)
  # 
  # #Find nearest positive definite matrix
  # corY = nearPD(corY,corr=TRUE,conv.tol=100000)$mat
  
  #corY = cor(Y)
  #browser()
  return(list(Y=Y,corY=corY))
}

sim_errors_and_bic = function(sim_func, fit1_func, fit2_func, num_sim){
  res = list()
  #browser()
  for(i in 1:num_sim){
    # if(i==64){
    #   browser()
    # }
    #browser()
    print(i)
    sim = sim_func(seed=i)
    res1 = fit1_func(sim)
    #browser()
    res2 = fit2_func(sim)
    #browser()
    res[[i]] = c(res1[[1]][1,], res1[[1]][2,], res1[[2]], res1[[4]],
                 res2[[1]][1,], res2[[1]][2,], res2[[2]], res2[[4]])
    #print(res[[i]])
  }
  #browser()
  res_df = as.data.frame(t(sapply(seq_along(res),function(s) res[[s]])))
  names(res_df)=c(paste("mae1",1:dim(res2[[1]])[2],sep="."),
               paste("rmse1",1:dim(res2[[1]])[2],sep="."),
               "bic1",paste0("param1.",1:length(res1[[4]])),
               paste("mae2",1:dim(res2[[1]])[2],sep="."),
               paste("rmse2",1:dim(res2[[1]])[2],sep="."),
               "bic2",paste0("param2.",1:length(res2[[4]])))
  return(res_df)
               
}
