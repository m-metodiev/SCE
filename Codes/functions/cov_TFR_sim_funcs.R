library(tmvtnorm)
sim_matList = function(dimension, beta, num_F, k_vec, num_G, F_0=TRUE){
  Fk = list()
  if(F_0){
    Fk[[1]] = diag(nrow=dimension,ncol=dimension)
  }
  for(k in (1:num_F)){
    X <- t(rmultinom(dimension, 1, rep(1/k_vec[k], k_vec[k]))) 
    Fk[[length(Fk) + 1]] = X%*%t(X)
  }
  
  Gl = list()
  Ml = list()
  Al = list()
  for(l in (1:num_G)){
    vois <- matrix(rbinom(dimension^2, 1, log(dimension)/dimension), dimension, dimension)
    vois <- (vois + t(vois)); diag(vois) <- 0; vois[which(vois > 1)] <- 1
    M = diag(rowSums(vois))
    A = vois
    A = (diag(1/rowSums(A))%*%A)
    Gl[[l]] = M - M%*%A*beta
    Ml[[l]] = M
    Al[[l]] = A
  }
  matList = list(Fk=Fk, Gl=Gl, Ml=Ml, Al=Al,vois=vois)
  return(matList)
}

# setting a seed for the bigger simulations
sim_cov_with_seed =function(num_observations, Sigma,seed=NULL,
                            degrees_of_freedom=NULL){

  dimension = dim(Sigma)[1]
  
  if(!is.null(seed)){
    set.seed(seed)
  }
  if(is.null(degrees_of_freedom)){
    dataset <- mvtnorm::rmvnorm(num_observations,
                                sigma=Sigma) + 2.1
  } else{
    #browser()
    # simulate from the t-distribution
    # multiply by a factor (does not change the estimation procedure)
    dataset <- tmvtnorm::rtmvt(num_observations,
                               mean=rep(2.1,nrow(Sigma)),
                               sigma=Sigma,
                               df=degrees_of_freedom,
                               lower=rep(0,nrow(Sigma)))*
      sqrt((degrees_of_freedom-2)/degrees_of_freedom)
  }
  
  correlation_matrix = cor_from_standard_errors(dataset-2.1)
  return(list(dataset=dataset,correlation_matrix=correlation_matrix))
}

sim_cov =function(num_observations, Sigma,seed=NULL, degrees_of_freedom=NULL){
  
  dimension = dim(Sigma)[1]
  if(is.null(degrees_of_freedom)){
    dataset <- mvtnorm::rmvnorm(num_observations, sigma=Sigma) + 2.1
  } else{
    #browser()
    # simulate from the t-distribution
    # multiply by a factor (does not change the estimation procedure)
    dataset <- tmvtnorm::rtmvt(num_observations,
                               mean=rep(2.1,nrow(Sigma)),
                               sigma=Sigma,
                               df=degrees_of_freedom,
                               lower=rep(0,nrow(Sigma))) *
      sqrt((degrees_of_freedom-2)/degrees_of_freedom)
  }
  
  
  # the dataset is already standardized
  correlation_matrix = cor_from_standard_errors(dataset-2.1)
  
  return(list(dataset=dataset,correlation_matrix=correlation_matrix))
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