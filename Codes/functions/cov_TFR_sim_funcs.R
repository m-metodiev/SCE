sim_matList = function(n, rho, num_F, k_vec, num_G, F_0=TRUE){
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
  return(list(Y=Y,corY=corY))
}

sim_cov =function(p, Sigma,seed=NULL){
  
  n = dim(Sigma)[1]
  Y <- mvtnorm::rmvnorm(p, sigma=Sigma)
  
  # Y is already standardized
  corY = cor_from_standard_errors(Y)
  
  return(list(Y=Y,corY=corY))
}

sim_errors_and_bic = function(sim_func, fit1_func, fit2_func, num_sim){
  res = list()
  
  for(i in 1:num_sim){
    print(i)
    sim = sim_func(seed=i)
    res1 = fit1_func(sim)
    res2 = fit2_func(sim)
    res[[i]] = c(res1[[1]][1,], res1[[1]][2,], res1[[2]], res1[[5]],
                 res2[[1]][1,], res2[[1]][2,], res2[[2]], res2[[5]])
  }
  
  res_df = as.data.frame(t(sapply(seq_along(res),function(s) res[[s]])))
  names(res_df)=c(paste("mae1",1:dim(res2[[1]])[2],sep="."),
               paste("rmse1",1:dim(res2[[1]])[2],sep="."),
               "bic1",paste0("param1.",1:length(res1[[5]])),
               paste("mae2",1:dim(res2[[1]])[2],sep="."),
               paste("rmse2",1:dim(res2[[1]])[2],sep="."),
               "bic2",paste0("param2.",1:length(res2[[5]])))
  return(res_df)
               
}