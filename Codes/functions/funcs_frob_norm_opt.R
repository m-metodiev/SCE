## Finding parameters that optimize the Frobenius norm

frob_scalar_prod = function(A,B){
  sum(diag(t(A)%*%B))
}

calc_Dmat = function(matList_full){
  Dmat = matrix(nrow=length(matList_full),ncol=length(matList_full))
  for(i in seq_along(matList_full)){
    for(j in seq_along(matList_full)){
      Dmat[i,j] = frob_scalar_prod(matList_full[[i]],matList_full[[j]])
    }
  }
  return(2*Dmat)
}

calc_dvec = function(S, matList_full){
  #browser()
  dvec = numeric(length(matList_full))
  for(i in seq_along(matList_full)){
    dvec[i] = frob_scalar_prod(S,matList_full[[i]])
  }
  return(2*dvec)
}

#a measure of distance between the support of two matrices
mat_support_distance = function(mat1, mat2){
  return(mean(abs(ceil(mat1)-ceil(mat2))))
}

#minimize the Frobenius norm via quadratic optimization
calc_Sigma_opt_frob = function(matList_full, covY, edge_constraints=c()){
  #browser()
  if(length(matList_full)==0){
    browser()
  }
  n = dim(matList_full[[1]])[1]
  #browser()
  #true_var_length = length(matList_full)
  
  # null_matrix_vec = rep(TRUE,length(matList_full))
  # for(i in (1:length(matList_full))){
  #   diag(matList_full[[i]]) = 0
  #   null_matrix_vec[i] = !(round(sum(matList_full[[i]]),15)==0)
  # }
  # matList_full = matList_full[null_matrix_vec]
  diag(covY)=0
  for(i in (1:length(matList_full))){
    diag(matList_full[[i]]) = 0
  } # we are only comparing the covariates, not the diagonals
  
  #The Frobenius inner product is included in each element
  Dmat = calc_Dmat(matList_full)
  dvec =  calc_dvec(covY,matList_full)
  num_param = length(dvec)
  
  #browser()
  #add constraint that sum has to be smaller than 1
  Dmat = rbind(cbind(Dmat,rep(0,dim(Dmat)[1])),c(rep(0,dim(Dmat)[1]),1))
  dvec = c(dvec,0)
  Amat = rbind(cbind(rep(-1,length(matList_full)),
                     diag(1,length(matList_full))),
               c(rep(0,length(matList_full)+1)))
  bvec = -t(t(c(1,numeric(length(matList_full)))))
  
  num_extra_var = 1 # some variables contribute nothing and only add constraints
  Sigma_0_opt = solve.QP(Dmat,dvec,Amat,bvec)
  init <- c(log(abs(Sigma_0_opt$solution)))[1:num_param]
  #browser()
  
  # need to stop parameters from being on the edge (messes up init)
  if(length(edge_constraints)>0){
    for(edge_constraint in edge_constraints){
      #browser()
      num_extra_var = num_extra_var + 1
      Dmat = rbind(cbind(Dmat,rep(0,dim(Dmat)[1])),c(rep(0,dim(Dmat)[1]),1))
      dvec = c(dvec,0)
      #browser()
      #pos_matrices = which(0!=round(sapply(seq_along(matList_full),function(s) exp(init)[s]*min(matList_full[[s]][matList_full[[s]]>0])),10))
      if(edge_constraint$r_min==0){
        #browser()
        #choose matrix closest to the identity matrix (with diagonal set to 0)
        #s_min = which.min(sapply(seq_along(matList_full),
        #                         function(s) mean(abs(ceiling(matList_full[[s]])-matrix(0,n,n))))[pos_matrices])
        #s_min = pos_matrices[s_min]
  
        constraint_vec = rep(-1,length(matList_full))
        #constraint_vec[s_min] = -2 # 1-sum(init) should be larger than init[s_min]
        constraint_digit = -1+edge_constraint$constraint_digit#exp(init)[edge_constraint$s_min]*min(matList_full[[edge_constraint$s_min]][matList_full[[edge_constraint$s_min]]>0])
      } else{
        #browser()
        #choose matrix closest to r_min matrix
        #r_min = edge_constraint
        #if(length(which(pos_matrices==r_min))>0){
        #  pos_matrices = pos_matrices[-which(pos_matrices==r_min)]
        #}
        #index = c(0,(1:num_param)[-r_min])
        #s_min = which.min(c(mean(abs(ceiling(matList_full[[r_min]])-matrix(0,n,n))),
        #                    sapply(seq_along(matList_full)[pos_matrices],
        #                           function(s) mean(abs(ceiling(matList_full[[s]])-
                                                          #ceiling(matList_full[[r_min]]))))))-1
        if(edge_constraint$s_min==0){
          #browser()
          constraint_vec = rep(0,length(matList_full))
          constraint_vec[edge_constraint$r_min] = 1
          constraint_digit = edge_constraint$constraint_digit
        } else{
          #browser()
          #s_min = pos_matrices[s_min]
          constraint_vec = rep(0,length(matList_full))
          constraint_vec[edge_constraint$r_min] = 1
          constraint_digit = edge_constraint$constraint_digit#exp(init)[edge_constraint$s_min]*min(matList_full[[edge_constraint$s_min]][matList_full[[edge_constraint$s_min]]>0])
        }
      }
      Amat = rbind(cbind(c(constraint_vec,
                           rep(0,dim(Amat)[1]-length(matList_full))),
                         Amat),
                   c(rep(0,dim(Amat)[1]+1)))
      bvec = t(t(c(constraint_digit,bvec)))
      Sigma_0_opt = solve.QP(Dmat,dvec,Amat,bvec)
      init <- c(log(abs(Sigma_0_opt$solution)))[1:num_param]
    }
  }
  value = Sigma_0_opt$value
  # init <- numeric(length(null_matrix_vec)) 
  # init[which(null_matrix_vec)] = res

  #browser()
  
  
  return(list(init=init,value=Sigma_0_opt$value))
}
