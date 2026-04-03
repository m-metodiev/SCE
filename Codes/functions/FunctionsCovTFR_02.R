library(gear)
library(ohenery)
library(tidyverse)
library(tictoc)
library(missMDA)

tilde_G_inv = function(M=NULL, A=NULL, beta,
                       U_full=NULL, solve_U_full=NULL,
                       solve_M_no_islands=NULL, eigen_real=NULL,
                       return_U_D_M=FALSE){
  # round down if too close to the edge
  # browser()
  if(beta==1){
    beta=.999
  }
  beta = beta*(beta<1-1e-4)+1e-4*(beta>=1-1e-4)
  
  # set diagonal to 1 for isolated nodes (A can't be defined for those)
  no_islands_id = diag(M)!=0
  G = matrix(0,ncol=dim(M)[1],nrow=dim(M)[1])
  
  # it is possible to save U_full to quicken the computation
  M = as.matrix(M)
  A = as.matrix(A)
  
  M_no_islands = M[no_islands_id,no_islands_id]
  A_no_islands = A[no_islands_id,no_islands_id]
  dimension = dim(A_no_islands)[1]
  
  # Matrix is only positive semidefinite (the are eigenvalues equal to 0),
  # due to the graph not being connected
  rankA = rankMatrix(A_no_islands)[1]
  
  if(is.null(U_full)){
    U_real = Re(eigen(A_no_islands)$vectors[,1:rankA])
    N_real = nullspace(A_no_islands)
    U_full = cbind(U_real,N_real)
    eigen_real = Re(eigen(A_no_islands)$values[1:rankA])
    
    D_beta_full = diag(c(eigen_real * beta / (1 - eigen_real * beta),rep(0,dimension-rankA)))
    
    G_inv_no_islands = (diag(dimension) + U_full%*%D_beta_full%*%solve(U_full))%*%solve(M_no_islands)
  } else{
    D_beta_full = diag(c(eigen_real * beta / (1 - eigen_real * beta),rep(0,dimension-rankA)))
    
    G_inv_no_islands = (diag(dimension) + U_full%*%D_beta_full%*%solve_U_full)%*%solve_M_no_islands
    
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

#' Computes correlation matrix from a normalized dataset (=standard errors)
#'
#' @param varepsilon n x d matrix (n = number of observations, d = dimension)
#'
#' @returns the correlation matrix
#' @keywords internal
cor_from_standard_errors = function(varepsilon){
  dimension=dim(varepsilon)[2]
  num_observations = dim(varepsilon)[1]
  correlation_matrix = matrix(0,ncol=dimension,nrow=dimension)
  for(t in 1:num_observations){
    correlation_matrix = correlation_matrix +
      t(t(varepsilon[t,]))%*%t(varepsilon[t,])
  }
  correlation_matrix = correlation_matrix/(num_observations-1)
  diag(correlation_matrix)=1
  library(Matrix)
  
  #Find nearest positive definite matrix
  correlation_matrix = nearPD(correlation_matrix,corr=TRUE,conv.tol=100000)$mat
  return(as.matrix(correlation_matrix))
}


CovMat_03 <- function(parm, matList, adj_positions, combined_effects=FALSE,
                      interaction_effects=list()){
  #browser()
  dimension <- ncol(matList$Fk[[1]])
  alpha <- parm[1:length(matList$Fk)]
  delta <- parm[length(matList$Fk)+(1:length(matList$Ml))]
  Phi <- Gamma <- matrix(0, dimension, dimension)
  
  if(combined_effects=="FosdickRaftery"){
    for(k in 1:length(matList$Fk)){Phi <- Phi + alpha[k]*matList$Fk[[k]]}
    adj_matrix = matList$Al[[1]]
    adj_matrix[is.na(adj_matrix)]=0
    adj_matrix = (adj_matrix != 0) + 0
    
    Gamma <- delta[1] * adj_matrix[adj_positions,adj_positions]
    Sigma <- .5*(Phi + Gamma + t(Phi+Gamma))#Phi + Gamma
    diag(Sigma) = 1
    return(list(Phi=Phi, Gamma=Gamma, Sigma=Sigma))
  }
  
  betal = parm[length(matList$Fk)+length(matList$Ml)+(1:length(matList$Ml))]
  s = dim(matList$Ml[[1]])[1]
  
  matList$Gl[[1]] = tilde_G_inv(matList$Ml[[1]],matList$Al[[1]],parm[length(parm)])[adj_positions,adj_positions]
  matList_combined = combined_matList(matList,interaction_effects=interaction_effects)$matList_full
  
  alpha_delta = parm[1:length(matList_combined)]
  Phi <- Gamma <- matrix(0, dimension, dimension)
  for(k in 1:length(matList_combined)){Phi <- Phi + alpha_delta[k]*matList_combined[[k]]}
  Sigma <- .5*(Phi + t(Phi));#Sigma is just a convex combination of all matrices
  
  diag(Sigma) = 1
  return(list(Phi=Phi, Gamma=Gamma, Sigma=Sigma, matList_combined=matList_combined,alpha_delta=alpha_delta))
}

CovMat_03 <- function(parm, matList, adj_positions, combined_effects=FALSE,
                      interaction_effects=list()){
  #browser()
  dimension <- ncol(matList$Fk[[1]])
  alpha <- parm[1:length(matList$Fk)]
  delta <- parm[length(matList$Fk)+(1:length(matList$Ml))]
  Phi <- Gamma <- matrix(0, dimension, dimension)
  
  if(combined_effects=="FosdickRaftery"){
    for(k in 1:length(matList$Fk)){Phi <- Phi + alpha[k]*matList$Fk[[k]]}
    adj_matrix = matList$Al[[1]]
    adj_matrix[is.na(adj_matrix)]=0
    adj_matrix = (adj_matrix != 0) + 0
    
    Gamma <- delta[1] * adj_matrix[adj_positions,adj_positions]
    Sigma <- .5*(Phi + Gamma + t(Phi+Gamma))#Phi + Gamma
    diag(Sigma) = 1
    return(list(Phi=Phi, Gamma=Gamma, Sigma=Sigma))
  }
  
  betal = parm[length(matList$Fk)+length(matList$Ml)+(1:length(matList$Ml))]
  s = dim(matList$Ml[[1]])[1]
  
  matList$Gl[[1]] = tilde_G_inv(matList$Ml[[1]],matList$Al[[1]],parm[length(parm)])[adj_positions,adj_positions]
  matList_combined = combined_matList(matList,interaction_effects=interaction_effects)$matList_full
  
  alpha_delta = parm[1:length(matList_combined)]
  Phi <- Gamma <- matrix(0, dimension, dimension)
  for(k in 1:length(matList_combined)){Phi <- Phi + alpha_delta[k]*matList_combined[[k]]}
  Sigma <- .5*(Phi + t(Phi));#Sigma is just a convex combination of all matrices
  
  diag(Sigma) = 1
  return(list(Phi=Phi, Gamma=Gamma, Sigma=Sigma, matList_combined=matList_combined,alpha_delta=alpha_delta))
}

#' Adds combined effects to the matList via the Hadamard product
#'
#' @inheritParams combined_matList_partial_der
#'
#' @returns named list with all matList combinations and their positions
#' @keywords internal
combined_matList = function(matList,interaction_effects=NULL,check_redundancy=FALSE){
  
  matList_full = c(matList$Fk,matList$Gl)
  for(i in seq_along(matList_full)){
    matList_full[[i]] = as.matrix(matList_full[[i]])
  }
  matList_full_names = c(names(matList$Fk),"spatial")
  counter = length(matList_full)
  sequence = seq_along(matList_full)
  matrix_pairs = list()
  found_redundant_pairs = FALSE
  
  #calculate all possible Hadamard-products; Exclude global effect matrix
  for(i in sequence){
    for(j in sequence[-(1:i)]){
      #browser()
      matprod = matList_full[[i]] * matList_full[[j]]
      if(check_redundancy){
        not_redundant = (Matrix::rankMatrix(sapply(c(matList_full,list(matprod)),
                                                   function(x) c(x)))[1] - (length(matList_full) + 1)) == 0
        if(!not_redundant){
          found_redundant_pairs = TRUE
        }
      } else{
        not_redundant = TRUE
      }
      if(!is.null(interaction_effects)){
        # only loop over those
        in_left = matList_full_names[i] == sapply(interaction_effects,function(interaction) interaction[1])
        in_both = sum(matList_full_names[j] == sapply(interaction_effects[in_left],function(interaction) interaction[2])) >= 1
        in_right = matList_full_names[j] == sapply(interaction_effects,function(interaction) interaction[1])
        in_both = in_both | (sum(matList_full_names[j] == sapply(interaction_effects[in_right],function(interaction) interaction[2]))>=1)
        include_matrix_pair = not_redundant & in_both
      } else{
        include_matrix_pair = not_redundant
      }
      if(include_matrix_pair){
        counter = counter + 1
        matList_full[[counter]] = matprod
        matrix_pairs = c(matrix_pairs,list(c(i,j)))
      }
    }
  }
  
  # add names of the matrices
  base_effects = c(names(matList$Fk),"spatial")
  effects = base_effects
  for(matrix_pair in matrix_pairs){
    effects[length(effects)+1] = paste0(base_effects[matrix_pair[1]]," and ",base_effects[matrix_pair[2]])
  }
  names(matList_full) = effects
  
  if(found_redundant_pairs){
    outputline = "Redundant effect pairs have been found. The following pairs are not redundant: "
    print(paste0(outputline,effects),quote=FALSE)
    print("Please use these effect pairs as input and start over.")
    return(-1)
  } else{
    return(list(matList_full=matList_full,matrix_pairs=matrix_pairs))
  }
}

#' Computes the matrices needed for the spatial effect
#'
#' @param adj_matrix the adjacency matrix of the spatial effect
#'
#' @returns a list containing Ml and Al, which are used to compute the
#'          correlation matrix for the spatial effect
get_M_A = function(adj_matrix){
  M = diag(rowSums(adj_matrix))
  A = adj_matrix
  A = (diag(1/rowSums(A))%*%A)
  Ml = list(M)
  Al = list(A)
  return(list(Ml=Ml,Al=Al))
}

