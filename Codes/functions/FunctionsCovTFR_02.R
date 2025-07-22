library(gear)
library(ohenery)
library(tidyverse)
library(tictoc)
library(missMDA)

calc_mean_neighbor_effect = function(matList, rho, beta, id_min){
  G_inv = calc_tilde_G_inv(matList$Ml[[1]],matList$Al[[1]],rho)[id_min,id_min]
  A = matList$Al[[1]][id_min,id_min]
  A[is.na(A)] = 0
  return(mean(beta*G_inv[as.matrix(A!=0)]))
}

forward_transform_param = function(param){
  
  # "Round down" if parameters are too close to the edge
  param[1:(length(param)-1)][param[1:(length(param)-1)]>=(1-1e-8)]=.99
  param[1:(length(param)-1)][param[1:(length(param)-1)]<1e-8]=rep(0.001/(length(param)-1),sum(param[1:(length(param)-1)]<1e-8))
  if(sum(param[1:(length(param)-1)])>=(1-1e-8)){
    param[1:(length(param)-1)] = param[1:(length(param)-1)]/sum(param[1:(length(param)-1)])*(1-1e-8)
  }
  param[length(param)] = param[length(param)]*(param[length(param)]<1-1e-4)+1e-4*(param[length(param)]>=1-1e-4)
  
  transformed_init = param
  transformed_init[1:(length(param)-1)] = inv_smax(c(1-sum(transformed_init[1:(length(param)-1)]),transformed_init[1:(length(param)-1)]))[-1]
  transformed_init[length(param)] = logit(transformed_init[length(param)])
  return(transformed_init)
}

backward_transform_param = function(param){
  transformed_parm = param
  transformed_parm[1:(length(param)-1)] = smax(c(-sum(transformed_parm[1:(length(param)-1)]),transformed_parm[1:(length(param)-1)]))[-1]
  transformed_parm[length(param)] = sigmoid(transformed_parm[length(param)])
  return(transformed_parm)
}

backward_transform_param_jacobian = function(param){
  
  # Jacobian of the softmax
  soft_max_der = function(par) diag(par)-
    matrix(rep(par,length(par)),ncol=length(par),nrow=length(par))*t(matrix(rep(par,length(par)),ncol=length(par),nrow=length(par)))
  
  # multivariate chain rule component
  sum_der = function(par) rbind(-1,diag(length(par)))
  
  jacobian_parm_full = param
  jacobian_parm = param[1:(length(param)-1)]
  jacobian_parm_matrix = cbind(0,diag(length(jacobian_parm)))%*%
    soft_max_der(smax(c(-sum(jacobian_parm),jacobian_parm)))%*%sum_der(jacobian_parm)

  jacobian_parm_full[length(param)] = sigmoid(jacobian_parm_full[length(param)])*(1-sigmoid(jacobian_parm_full[length(param)])) # derivative of the sigmoid function
  jacobian_parm_matrix = rbind(cbind(jacobian_parm_matrix,rep(0,length(param)-1)),c(rep(0,length(param)-1),jacobian_parm_full[length(param)]))
  return(jacobian_parm_matrix)
}

calc_tilde_G_inv = function(M=NULL, A=NULL, rho, 
                            U_full=NULL, solve_U_full=NULL, 
                            solve_M_no_islands=NULL, eigen_real=NULL,
                            return_U_D_M=FALSE){
  # round down if too close to the edge
  if(rho==1){
    rho=.999
  }
  rho = rho*(rho<1-1e-4)+1e-4*(rho>=1-1e-4)

  # set diagonal to 1 for isolated nodes (A can't be defined for those)  
  no_islands_id = diag(M)!=0
  G = matrix(0,ncol=dim(M)[1],nrow=dim(M)[1])
  
  # it is possible to save U_full to quicken the computation
  M = as.matrix(M)
  A = as.matrix(A)
  
  M_no_islands = M[no_islands_id,no_islands_id]
  A_no_islands = A[no_islands_id,no_islands_id]
  n = dim(A_no_islands)[1]
  
  # Matrix is only positive semidefinite (the are eigenvalues equal to 0), 
  # due to the graph not being connected
  rankA = rankMatrix(A_no_islands)[1]
  
  if(is.null(U_full)){
    U_real = Re(eigen(A_no_islands)$vectors[,1:rankA])
    N_real = nullspace(A_no_islands)
    U_full = cbind(U_real,N_real)
    eigen_real = Re(eigen(A_no_islands)$values[1:rankA])
    
    D_rho_full = diag(c(eigen_real * rho / (1 - eigen_real * rho),rep(0,n-rankA))) 
    
    G_inv_no_islands = (diag(n) + U_full%*%D_rho_full%*%solve(U_full))%*%solve(M_no_islands)
  } else{
    D_rho_full = diag(c(eigen_real * rho / (1 - eigen_real * rho),rep(0,n-rankA))) 

    G_inv_no_islands = (diag(n) + U_full%*%D_rho_full%*%solve_U_full)%*%solve_M_no_islands
    
  }

  G[no_islands_id,no_islands_id] = cov2cor(G_inv_no_islands)
  diag(G)=1 # islands are independent of all other nodes
  
  if(return_U_D_M){ # to compute the real value
    return(list(U_full=U_full, solve_U_full=solve(U_full), 
                solve_M_no_islands=solve(M_no_islands), eigen_real=eigen_real))
  } else{
    return(G)
  }
}

calc_tilde_G_inv_partial_rho = function(M, A, rho){
  
  M = as.matrix(M)
  A = as.matrix(A)
  
  # set diagonal to 1 for isolated nodes (A can't be defined for those)
  no_islands_id = diag(M)!=0
  G_inv_partial_rho = matrix(0,ncol=dim(M)[1],nrow=dim(M)[1])
  G_inv = diag(dim(M)[1])
  
  # Matrix is only positive semidefinite due to the graph not being connected
  M_no_islands = M[no_islands_id,no_islands_id]
  A_no_islands = A[no_islands_id,no_islands_id]
  n = dim(A_no_islands)[1]
  rankA = rankMatrix(A_no_islands)[1]
  U_real = Re(eigen(A_no_islands)$vectors[,1:rankA])
  N_real = nullspace(A_no_islands)
  U_full = cbind(U_real,N_real)
  eigen_real = Re(eigen(A_no_islands)$values[1:rankA])
  D_rho_full = diag(c(eigen_real * rho / (1 - eigen_real * rho),rep(0,n-rankA))) 
  D_rho_full_partial_rho = diag(c(eigen_real/ ((1 - eigen_real * rho)^2),rep(0,n-rankA)))
  solve_U = solve(U_full)
  solve_M = solve(M_no_islands)
  
  G_inv_partial_rho_no_islands = (U_full%*%D_rho_full_partial_rho%*%solve_U)%*%solve_M
  G_inv_partial_rho[no_islands_id,no_islands_id] = G_inv_partial_rho_no_islands
  G_inv_no_islands = (diag(n) + U_full%*%D_rho_full%*%solve_U)%*%solve_M
  G_inv[no_islands_id, no_islands_id] = G_inv_no_islands
  
  S_inv_sqrt_no_islands = diag(1/sqrt(diag(G_inv_no_islands)))
  S_inv_sqrt_partial_rho_no_islands = (-1/2)*diag(1/(diag(G_inv_no_islands)^(3/2))*diag(G_inv_partial_rho_no_islands))
  
  tilde_G_inv_partial_rho_no_islands = S_inv_sqrt_partial_rho_no_islands %*% G_inv_no_islands %*% S_inv_sqrt_no_islands+
    S_inv_sqrt_no_islands %*% G_inv_partial_rho_no_islands %*% S_inv_sqrt_no_islands +
    S_inv_sqrt_no_islands %*% G_inv_no_islands %*% S_inv_sqrt_partial_rho_no_islands
  
  #islands are constant, so their derivatives are all 0
  tilde_G_inv_partial_rho = matrix(0,ncol=dim(M)[1],nrow=dim(M)[1])
  tilde_G_inv_partial_rho[no_islands_id, no_islands_id] = tilde_G_inv_partial_rho_no_islands

  # the diagonal values are also constant, so the diagonal derivative is also 0
  diag(tilde_G_inv_partial_rho) = 0
  
  # return the correlation matrix at the end
  G_inv[no_islands_id, no_islands_id] = cov2cor(G_inv_no_islands)
  return(list(tilde_G_inv= G_inv,tilde_G_inv_partial_rho=tilde_G_inv_partial_rho))
}

# the other loglikelihood is translated; I'm just including the normalizing
# constant for the BIC;
true_LogLikParm_02 <- function(id_min, parm, matList, Y,
                          link=function(matList) c(matList$Fk,matList$Gl)){
  num_observations = nrow(Y)
  has_missing_values = sum(is.na(Y))>0
  n = dim(matList$Ml[[1]])[1]
  res = 0
  
  for(t in (1:num_observations)){
    id_min_t = !is.na(Y[t,])
    res = res + mvtnorm::dmvnorm(Y[t,id_min_t], sigma=as.matrix(CovMat_03(id_min=id_min, parm=parm, matList=matList, link=link)$Sigma)[id_min_t,id_min_t], log=TRUE)
  }
  return(res)
}

LogLikParm_02 <- function(id_min, parm, matList, Y,
                          link=function(matList) c(matList$Fk,matList$Gl)){
  #browser()
  num_observations = nrow(Y)
  has_missing_values = sum(is.na(Y))>0
  if(has_missing_values){
    n = dim(matList$Ml[[1]])[1]
    res = 0
    
    for(t in (1:num_observations)){
      #browser()
      id_min_t = !is.na(Y[t,])
      Sigma = as.matrix(CovMat_03(id_min=id_min, 
                                  parm=parm, matList=matList, 
                                  link=link)$Sigma)[id_min_t,id_min_t]
      S = Y[t,id_min_t]%*%t(Y[t,id_min_t])
      res = res -sum(log(eigen(Sigma)$values))-sum(diag(S%*%solve(Sigma)))
    }
    
  } else{
    Sigma = as.matrix(CovMat_03(id_min=id_min, 
                                parm=parm, matList=matList, 
                                link=link)$Sigma)
    S = Y[1,]%*%t(Y[1,])
    if(num_observations!=1){
      for(t in (2:num_observations)){
        S = S + Y[t,]%*%t(Y[t,])
      }
    }
    S = S/num_observations
    #browser()
    res = -sum(log(eigen(Sigma)$values))-sum(diag(S%*%solve(Sigma)))
  }
  return(res)
}

LogLikLogParm_02 <- function(id_min, logParm, matList, Y,
                             link=function(matList) c(matList$Fk,matList$Gl)){
  parm = backward_transform_param(logParm)
  # "Round down" if parameters are too close to the edge
  parm[1:(length(parm)-1)][parm[1:(length(parm)-1)]>=(1-1e-8)]=.99
  parm[1:(length(parm)-1)][parm[1:(length(parm)-1)]<1e-8]=rep(0.001/(length(parm)-1),sum(parm[1:(length(parm)-1)]<1e-8))
  if(sum(parm[1:(length(parm)-1)])>=(1-1e-8)){
    parm[1:(length(parm)-1)] = parm[1:(length(parm)-1)]/sum(parm[1:(length(parm)-1)])*(1-1e-8)
  }
  parm[length(parm)] = parm[length(parm)]*(parm[length(parm)]<1-1e-4)+1e-4*(parm[length(parm)]>=1-1e-4)
  LogLikParm_02(id_min=id_min, parm, matList, Y, link=link)
}

#adds combined effects to the matList via the Hadamard product
combined_matList = function(matList){
  
  matList_full = c(matList$Fk,matList$Gl)
  counter = length(matList_full)
  sequence = seq_along(matList_full)
  #calculate all possible Hadamard-products; Exclude global effect matrix
  for(i in sequence[-length(matList$Fk)]){
    for(j in sequence[c(-(1:i),-length(matList$Fk))]){
      counter = counter + 1
      matList_full[[counter]] = matList_full[[i]] * matList_full[[j]]
    }
  }
  return(matList_full)
}

#adds combined effects to the matList via the Hadamard product
link_der_combined = function(matList, link_matList, tilde_G_inv_partial_rho){
  #browser()
  n = dim(link_matList[[1]])[1]
  matList_full = c(matList$Fk,matList$Gl)
  counter = length(matList_full)
  sequence = seq_along(matList_full)
  #calculate all possible Hadamard-products; Exclude global effect matrix
  for(i in sequence){
    matList_full[[i]] = matrix(0,n,n) # does not depend on rho
  }
  
  matList_full[[length(matList$Fk)+1]] = tilde_G_inv_partial_rho
  
  for(i in sequence[-length(matList$Fk)]){
    for(j in sequence[c(-(1:i),-length(matList$Fk))]){
      
      counter = counter + 1
      
      # use product rule whenever tilde_G_inv is included
      if(i==(length(matList$Fk)+1)){
        matList_full[[counter]] = tilde_G_inv_partial_rho * link_matList[[j]]
      } else if(j==(length(matList$Fk)+1)){
        matList_full[[counter]] = link_matList[[i]] * tilde_G_inv_partial_rho
      } else{
        matList_full[[counter]] = matrix(0,n,n)
      }
    }
  }
  #browser()
  return(matList_full)
}

# list of derivatives for the matList wrt rho; could change depending on link
link_der_simple = function(matList, link_matList, tilde_G_inv_partial_rho){
  n = dim(link_matList[[1]])[1]
  matList_der = link_matList
  for(k in seq_along(link_matList)){
    matList_der[[k]] = matrix(0,n,n) #Fks do not depend on rho
  }
  matList_der[[length(matList_der)]] = tilde_G_inv_partial_rho
  return(matList_der)
}

# These are actually 2 completely different function:
# One returns the gradient of the loglikelihood, one the gradient of Sigma
# You can choose which one to use by setting the parameter return_Sigma_der
# Warning: If Y contains missing values, it is expected to be a matrix with only 1 row
GradLogLikParm_02 <- function(id_min, parm, matList, Y, 
                              link=function(matList) c(matList$Fk,matList$Gl),
                              link_der_rho=link_der_simple,
                              return_Sigma_der=FALSE){
  #browser()
  n = dim(matList$Fk[[1]])[1]
  p = nrow(Y)
  if(p==1){
    covY <- t(t(Y[,!is.na(Y[1,])]))%*%t(Y[,!is.na(Y[1,])])/p
  } else{
    covY <- t(Y)%*%Y/p
  }
  
  # Calculate derivatives for matrix from the CAR model
  l=1
  rhol = parm[length(parm)]
  G_inv_list = calc_tilde_G_inv_partial_rho(matList$Ml[[l]],matList$Al[[l]],rhol[l])
  G_inv_list$tilde_G_inv = G_inv_list$tilde_G_inv[id_min,id_min]
  G_inv_list$tilde_G_inv_partial_rho = G_inv_list$tilde_G_inv_partial_rho[id_min,id_min]
  matList$Gl[[1]] = G_inv_list$tilde_G_inv
  
  link_matList = link(matList) # all operations are performed in this matList
  gradLogLik.alpha_beta <- numeric(length(link_matList))
  
  gradLogLik.rho <- numeric(length(matList$Gl)) 
  # probably only works if numbers of Gl-matrices is equal to 1
  
  Sigma = as.matrix(CovMat_03(parm, matList,id_min=id_min,link=link)$Sigma)[!is.na(Y[1,]),!is.na(Y[1,])]
  
  # due to numeric issues, Sigma could be computationally singular
  # derivative is set to be very low if that happens
  if(sum(is.na(eigen(Sigma)$values))==0){
    Omega = solve(Sigma) 
  } else{
    Omega = -exp(16)*diag(n)
  }

  for(k in seq_along(link_matList)){
    diag(link_matList[[k]])=0 # there is an identity matrix as the first F
    gradLogLik.alpha_beta[k] = - sum(diag((link_matList[[k]][!is.na(Y[1,]),!is.na(Y[1,])])%*%Omega)) +
      sum(diag(covY%*%Omega%*%(link_matList[[k]][!is.na(Y[1,]),!is.na(Y[1,])])%*%Omega))
  }
  
  tilde_G_inv_partial_rho = G_inv_list$tilde_G_inv_partial_rho
  rho_der_list = link_der_rho(matList, link_matList, tilde_G_inv_partial_rho)
  Sigma_partial_rho = matrix(0,n,n)[!is.na(Y[1,]),!is.na(Y[1,])]
  for(k in seq_along(rho_der_list)){
    Sigma_partial_rho = Sigma_partial_rho + parm[k] * rho_der_list[[k]][!is.na(Y[1,]),!is.na(Y[1,])]
  }
  gradLogLik.rho[l] = (-sum(diag(Sigma_partial_rho%*%Omega))+sum(diag(covY%*%Omega%*%Sigma_partial_rho%*%Omega)))
  if(return_Sigma_der){
    return(c(link_matList=link_matList,list(Sigma_partial_rho=Sigma_partial_rho)))
  } else{
    return(c(gradLogLik.alpha_beta, gradLogLik.rho))
  }
}

Fisher_information = function(id_min, parm, matList, link, link_der_rho){
  #browser()
  Sigma = CovMat_03(parm, matList,id_min=id_min,link=link)$Sigma
  Sigma_inv = solve(Sigma)
  Sigma_der = GradLogLikParm_02(id_min, parm, matList, Y=matrix(0,ncol=dim(Sigma)[1],nrow=dim(Sigma)[1]), link=link,
                                link_der_rho=link_der_rho, return_Sigma_der=TRUE)
  Fisher_mat = matrix(0,ncol=length(parm),nrow=length(parm))

  for(i in seq_along(parm)){
    for(j in seq_along(parm)){
      Fisher_mat[i,j] = (1/2)*sum(diag(Sigma_inv%*%Sigma_der[[i]]%*%Sigma_inv%*%Sigma_der[[j]]))
    }
  }
  return(Fisher_mat) 
}

cor_from_standard_errors = function(varepsilon){
  n=dim(varepsilon)[2]
  p=dim(varepsilon)[1]
  corY = matrix(0,ncol=n,nrow=n)
  for(t in 1:p){
    corY = corY + t(t(varepsilon[t,]))%*%t(varepsilon[t,])
  }
  corY = corY/(p-1)
  diag(corY)=1
  library(Matrix)
  
  #Find nearest positive definite matrix
  corY = nearPD(corY,corr=TRUE,conv.tol=100000)$mat
  return(as.matrix(corY))
}

compute_marginal_cor = function(Y){
  n = ncol(Y)
  corY = matrix(ncol=n,nrow=n)
  
  #use pairwise correlation estimates
  for(i in (1:n)){
    for(j in (i:n)){
      #browser()
      Yi_notmissing = !is.na(Y[,i])
      Yj_notmissing = !is.na(Y[,j])
      corY_ij = cor(Y[Yi_notmissing&Yj_notmissing,i], Y[Yi_notmissing&Yj_notmissing,j])
      corY[i,j] = corY_ij
      corY[j,i] = corY_ij
    }
  }
  return(corY)
}

eta_D_der = function(parm, matList, id_min, link, link_der_rho,index=4){
  #browser()
  parm=sapply(c(parm),function(p)p)
  covMatstuff = CovMat_03(parm=parm ,matList=matList,id_min=id_min,link=link)
  ml_combined = covMatstuff$matList_combined
  Sigma = CovMat_03(parm, matList,id_min=id_min,link=link)$Sigma
  Sigma_der = GradLogLikParm_02(id_min, parm, matList, Y=matrix(0,ncol=dim(Sigma)[1],nrow=dim(Sigma)[1]), link=link,
                                link_der_rho=link_der_rho, return_Sigma_der=TRUE)
  
  # determine support
  matList_supp = matList
  
  #Set to 1 if countries are not neighbors
  matList_supp$Al[[1]][is.na(matList_supp$Al[[1]])] = 0
  matList_supp$Gl[[1]] = (matList_supp$Al[[1]][id_min,id_min] != 0) + 0
  ml_combined_supp = link(matList_supp)
  
  #diagonals are not in the sum
  diag(ml_combined[[index]]) = 0
  diag(ml_combined_supp[[index]]) = 0
  
  eta_D_der = c(sum(parm[index]*ml_combined[[index]]*ml_combined_supp[[index]]),
                sum(Sigma_der[[index]]*ml_combined_supp[[index]]))/sum(ml_combined_supp[[index]])
  return(eta_D_der)
}

#Calculate average effect (=mean effect over matrix support)
avg_effect = function(parm, matList, id_min, link){
  
  #Determine support of the matrices
  matList_supp = matList
  
  #Set to 1 if countries are not neighbors
  matList_supp$Al[[1]][is.na(matList_supp$Al[[1]])] = 0
  matList_supp$Gl[[1]] = (matList_supp$Al[[1]][id_min,id_min] != 0) + 0
  ml_combined_supp = link(matList_supp) 
  
  covMatstuff = CovMat_03(parm=parm ,matList=matList,id_min=id_min,link=link)
  ml_combined = covMatstuff$matList_combined
  alpha_beta = covMatstuff$alpha_beta
  
  for(Fk in ml_combined_supp){
    diag(Fk)=0
  }
  for(Fk in ml_combined){
    diag(Fk)=0
  }
  # comparing the diagonals makes no sense
  
  #browser()
  sapply(seq_along(alpha_beta), function(i) alpha_beta[i]*sum(ml_combined_supp[[i]]*ml_combined[[i]])/sum(ml_combined_supp[[i]]))
}

compute_lambda_opt = function(id_min, parm, matList, link, link_der_rho, Y, pearson_mat, SCE_mat,
                              Y_nonmissing=NULL, use_bootstrap=FALSE){
  
  # numerical issues
  if(parm[length(parm)]<1e-1){
    parm[length(parm)]=1e-1
  }
  if(parm[length(parm)-1]<1e-4){
    parm[length(parm)-1]=1e-4
  }
  
  # the covariance matrix is given by the inverse of the Fisher information
  Fisher_mat = Fisher_information(id_min, parm, matList, link, link_der_rho)
  Fisher_mat = solve(Fisher_mat)
  
  n = dim(matList$Fk[[1]])[1]
  p = dim(Y)[1]
  Sigma_der = GradLogLikParm_02(id_min, parm, matList, Y=matrix(0,ncol=n,nrow=n), link=link,
                                link_der_rho=link_der_rho, return_Sigma_der=TRUE)
  n_small = n
  total_der = matrix(0, ncol=n_small*n_small,nrow=length(parm))
  for(i in seq_along(parm)){
    total_der[i,] = c(as.matrix(Sigma_der[[i]]))
  }
  
  #compute approximate covariance matrix of the SCE by the delta method
  approx_var_SCE = sum(diag((t(total_der)%*%Fisher_mat%*%total_der)))
  
  # more expensive, but gives more accurate estimates
  if(sum(is.na(Y))>0 | use_bootstrap){
    test_func = function(s){
      set.seed(s) # ensures replicability
      test_data = mvtnorm::rmvnorm(p, sigma=pearson_mat)
      test_data[is.na(Y)]=NA
      
      # impute data if values are missing, use Pearson matrix otherwise
      if(sum(is.na(Y))>0){
        sim_test = list(Y=test_data,corY=compute_marginal_cor(test_data))
        pearson_test = cor_from_standard_errors(imputePCA(test_data,ncp=length(parm))$completeObs)#cor(imputePCA(test_data,ncp=length(parm))$completeObs)
      } else{
        sim_test = list(Y=test_data, corY=cor_from_standard_errors(test_data))
        pearson_test = cor(test_data)
      }
      SCE_test = fit_param(n, p, matList, sim_test,id_min=id_min, 
                           link=link, link_der_rho=link_der_rho, save=FALSE,
                           init=log(parm), simple=TRUE)
      SCE_test = SCE_test[[2]][[length(SCE_test[[2]])]]
      
      system(sprintf('echo "\n%s\n"', paste0(s, collapse="")))
      return(array( c(SCE_test , pearson_test ) , dim = c( n , n , 2 ) ))
    }
    #browser()
    pearson_mat=as.matrix(pearson_mat)
    num_bootstrap_iters = 100
    bootstrap_sample = mclapply(1:num_bootstrap_iters, FUN=function(s) test_func(s), mc.cores=5)
    bootstrap_sample = simplify2array(bootstrap_sample)

    (test_var = sum(apply(bootstrap_sample,1, function(bootstrap_sample_row) apply(bootstrap_sample_row,1, function(bootstrap_sample_row_col) var(bootstrap_sample_row_col[2,])*p)))) 
    (test_cov = sum(apply(bootstrap_sample,1, function(bootstrap_sample_row) apply(bootstrap_sample_row,1, function(bootstrap_sample_row_col) cov(bootstrap_sample_row_col[1,],bootstrap_sample_row_col[2,])*p))))
    (test_mse = sum((apply(bootstrap_sample,1, function(bootstrap_sample_row) apply(bootstrap_sample_row,1, function(bootstrap_sample_row_col) mean(bootstrap_sample_row_col[1,])))-pearson_mat)^2)) 
    (test_var2 = sum(apply(bootstrap_sample,1, function(bootstrap_sample_row) apply(bootstrap_sample_row,1, function(bootstrap_sample_row_col) var(bootstrap_sample_row_col[1,])*p))))   

    approx_pi =   test_mse
    approx_mse = sum((SCE_mat- pearson_mat )^2)
    (res = (test_var-test_cov)/(p*test_mse))
  } else{
    approx_pi = sum((1-pearson_mat^2)^2)
    approx_mse = sum((SCE_mat- pearson_mat )^2)
    approx_var_SCE = sum(diag(t(total_der)%*%Fisher_mat%*%total_der))
    (res = (sqrt(approx_pi)*(sqrt(approx_pi)-sqrt(approx_var_SCE)))/(p*approx_mse))
  }
  res = min(res,1)
  res = max(res,0)
  #lambda has to be between 0 and 1
  
  return(res)
}

GradLogLikLogParm_02 <- function(id_min, logParm, matList, Y, 
                                 link=function(matList) c(matList$Fk,matList$Gl),
                                 link_der_rho=link_der_simple){
  parm = backward_transform_param(logParm)
  
  # round down if parameter is on the edge
  parm[1:(length(parm)-1)][parm[1:(length(parm)-1)]>=(1-1e-8)]=.99
  parm[1:(length(parm)-1)][parm[1:(length(parm)-1)]<1e-8]=rep(0.001/(length(parm)-1),sum(parm[1:(length(parm)-1)]<1e-8))
  parm[length(parm)] = parm[length(parm)]*(parm[length(parm)]<1-1e-4)+1e-4*(parm[length(parm)]>=1-1e-4)
  
  jacobian = backward_transform_param_jacobian(logParm)
  
  if(sum(is.na(Y))==0){
    gradient=GradLogLikParm_02(id_min=id_min, parm=parm, matList, Y, link=link,link_der_rho=link_der_rho)
  } else{
    # use linearity of the determinant
    gradient=rowSums(sapply(1:nrow(Y),function(t) GradLogLikParm_02(id_min=id_min, parm=parm, matList, t(Y[t,]),
                                                                    link=link,link_der_rho=link_der_rho)))
  }
  return(gradient%*%jacobian)
}



CovMat_03 <- function(parm, matList, id_min, combined_effects=FALSE,
                      link=function(matList) c(matList$Fk,matList$Gl)){
  
  n <- ncol(matList$Fk[[1]])
  alpha <- parm[1:length(matList$Fk)]
  beta <- parm[length(matList$Fk)+(1:length(matList$Ml))]
  Phi <- Gamma <- matrix(0, n, n)
  
  if(combined_effects=="FosdickRaftery"){
    for(k in 1:length(matList$Fk)){Phi <- Phi + alpha[k]*matList$Fk[[k]]}
    adj_matrix = matList$Al[[1]]
    adj_matrix[is.na(adj_matrix)]=0
    adj_matrix = (adj_matrix != 0) + 0
    
    Gamma <- beta[1] * adj_matrix[id_min,id_min]
    Sigma <- .5*(Phi + Gamma + t(Phi+Gamma))#Phi + Gamma
    diag(Sigma) = 1
    return(list(Phi=Phi, Gamma=Gamma, Sigma=Sigma))
  }

  rhol = parm[length(matList$Fk)+length(matList$Ml)+(1:length(matList$Ml))]
  s = dim(matList$Ml[[1]])[1]
  
  matList$Gl[[1]] = calc_tilde_G_inv(matList$Ml[[1]],matList$Al[[1]],parm[length(parm)])[id_min,id_min]
  matList_combined = link(matList)
  
  alpha_beta = parm[1:length(matList_combined)]
  Phi <- Gamma <- matrix(0, n, n)
  for(k in 1:length(matList_combined)){Phi <- Phi + alpha_beta[k]*matList_combined[[k]]}
  Sigma <- .5*(Phi + t(Phi));#Sigma is just a convex combination of all matrices

  diag(Sigma) = 1
  return(list(Phi=Phi, Gamma=Gamma, Sigma=Sigma, matList_combined=matList_combined,alpha_beta=alpha_beta))
}