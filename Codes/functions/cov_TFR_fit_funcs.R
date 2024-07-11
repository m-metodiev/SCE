PACKAGES_MISC = c("sna", "xts") # stuff
PACKAGES_STAT_OPTIM = c("mvtnorm","quadprog","pracma", "Matrix","missMDA", "covFactorModel",
                        "quadprog", "cvCovEst", "CovTools", "CVglasso") # stat,math stuff
#PACKAGES_VIS = c("UpSetR", "ggplot2", "grid", "plyr", "igraph") # plot stuff
PACKAGES_PARALLEL = c("foreach", 
                      "parallel", "doParallel") # parallelization stuff
PACKAGES = c(PACKAGES_MISC, PACKAGES_STAT_OPTIM, 
             #PACKAGES_VIS, 
             PACKAGES_PARALLEL) # all of the stuff
lapply(PACKAGES, require, character.only = TRUE)
# library(sna); library(mvtnorm); library(quadprog); 
# library(igraph)
# library(cvCovEst); library(CovTools)
# 
# # for upset plot
# library(UpSetR)
# library(ggplot2)
# library(grid)
# library(plyr)
# 
# # for parallelization
# library(foreach)
# library(parallel)
# library(doParallel)

fit_param = function(n, p, matList, sim, id_min,
                     filename_param_fit=NULL, filename_ests=NULL, filename_bic=NULL, 
                     filename_error_measures=NULL,
                     num_grid_sim = 100, ncores=8,
                     save=TRUE,Sigma=NULL,
                     link=function(matList) c(matList$Fk,matList$Gl), 
                     link_der_rho=link_der_simple,
                     normalize_data=FALSE, compute_WSCE=FALSE,
                     simple=FALSE, init=NULL){
  #browser()
  
  matList$Gl[[1]] = matrix(0,n,n) # Just a dummy-variable, will be declared later
  p_is_one=FALSE
  if(p==1){
    sim$Y = t(matrix(sim$Y))
    p_is_one=TRUE
  } # code is slightly different if p==1
  #browser()
  # the SCE needs to be computed on the normalized dataset
  if(normalize_data){
    #browser()
    simY_normalized = (sim$Y - t(matrix(rep(colMeans(sim$Y),p),ncol=p)))%*%diag(1/apply(sim$Y,2,sd))
    sim$corY=cor(sim$Y)
    #simY_normalized = (sim$Y - t(matrix(rep(rep(mean(sim$Y),n),p),ncol=p)))%*%diag(1/sd(sim$Y))
  } else{
    simY_normalized = sim$Y
  }
  # for the real dataset, we know that the standardized errors are 
  # already normalized, so we don't need to do that again
  if(is.null(init)){
    init = fit_init(id_min,matList = matList, corY=sim$corY,
                    num_grid_sim=num_grid_sim, ncores=ncores,
                    link=link)
  }

  #browser()
  exp_param_optim_frob = exp(init)
  Sigma_0_optim_frob = CovMat_03(parm=exp_param_optim_frob, matList,id_min=id_min)$Sigma
  #browser()
  list2env(fit_cov(id_min, p,
                      matList,
                      forward_transform_param(exp(init)),
                      as.matrix(simY_normalized),
                        link=link,link_der_rho=link_der_rho), envir = .GlobalEnv)
  #browser()
  if(sum(is.na(sim$Y))>0){
    has_missingvalues=TRUE
    Y_nonmissing = imputePCA(sim$Y,ncp=length(init))$completeObs
  } else{
    has_missingvalues=FALSE
    Y_nonmissing = sim$Y
  }
  
  if((p>1) & (!simple)){
    #browser()
    ests = list(as.matrix(cor(Y_nonmissing)),
                cov2cor(covFactorModel(xts(x=Y_nonmissing, order.by=Sys.Date()-1:p),K=length(init))),
                cov2cor(CVglasso(X=Y_nonmissing)$Sigma),#cov2cor(solve(PreEst.glasso(sim$Y)$C)),
                cov2cor(linearShrinkLWEst(Y_nonmissing)),
                as.matrix(Sigma_0_optim_frob),
                as.matrix(SigmaHat1))
    
    if(compute_WSCE){
      #browser()
      if(normalize_data){
        pearson_mat = cor(sim$Y) # mean and variance unknown
        use_bootstrap = TRUE
      } else{
        pearson_mat = cor_from_standard_errors(Y_nonmissing) # mean and variance known
        use_bootstrap = FALSE # used anyway if part of the data is missing
      }
      #browser()
      lambda = compute_lambda_opt(id_min,parm=param_fit1,matList=matList,
                                  link=link,link_der_rho=link_der_rho,Y=sim$Y,
                                  pearson_mat=pearson_mat,SCE_mat=as.matrix(SigmaHat1),
                                  Y_nonmissing=Y_nonmissing, use_bootstrap=use_bootstrap)
      #browser()
      print("lambda")
      print(lambda)
      ests = c(ests, list(lambda*as.matrix(SigmaHat1)+(1-lambda)*as.matrix(cor(Y_nonmissing))))
    }
  } else{
    #browser()
    ests = list(as.matrix(Sigma_0_optim_frob),as.matrix(SigmaHat1))
  } # If p=1, most estimators cant be computed  

  #browser()
  
  if(save==TRUE){
    #browser()
    write_param(t(c(param_fit1)),filename=filename_param_fit,
                matList = matList, link=link)
    write.csv(sapply(seq_along(ests),function(s) ests[[s]]),file=filename_ests)
    write.csv(as.data.frame(bic),file=filename_bic)
    if(!is.null(filename_error_measures)){
      #browser()
      write_summary_measures(filename_ests = filename_ests, 
                             filename_error_measures = filename_error_measures,
                             p_is_one=p_is_one,Sigma=Sigma)
    }
  } else{
    if(!is.null(filename_error_measures)){
      summary_measures = sapply(1:length(ests),
                                function(i) c(mean(abs(ests[[i]]-Sigma)),
                                              sqrt(mean((ests[[i]]-Sigma)^2))))
      return(list(summary_measures, bic, ests, param_fit1))
    } else{ # error measures aren't necessarily always available
      return(list(t(c(param_fit1)), ests, bic))
    }
  }
  #browser()
  
}

fit_init = function(id_min, matList, corY,
                    num_grid_sim = 100, ncores=8,
                    link=function(matList) c(matList$Fk,matList$Gl)){
  #browser()
  n = dim(matList$Fk[[1]])[1]
  s = dim(matList$Al[[1]])[1]
  #Need grid-search for rho because the norm is not quadratic w.r.t. rho
  beta=1.5
  xi = (1:(num_grid_sim+1))/(num_grid_sim+1)
  #tan-hyperbolic-spaced grid because rho approaches 1
  rho_vec = (1-tanh(beta*(1+xi))/tanh(beta))/(min((1-tanh(beta*(1+xi))/tanh(beta))))
  rho_vec = rho_vec[-length(rho_vec)]
  
  is_on_edge = TRUE # solution can lie on the edge of the parameter space
  edge_constraints = list() # params which lie on the edge will be adjusted in constraints
  
  counter = 0
  while(is_on_edge){
    #browser()
    counter = counter + 1
    grid_search = function(rho){
      #browser()
      matList$Gl[[1]] = calc_tilde_G_inv(M=matList$Ml[[1]],A=matList$Al[[1]],rho=rho,
                                         U_full=matList$U_full, solve_U_full=matList$solve_U_full,
                                         solve_M_no_islands=matList$solve_M_no_islands,
                                         eigen_real=matList$eigen_real)[id_min,id_min]
      matList_full = link(matList)
      #browser()
      res = calc_Sigma_opt_frob(matList_full, corY, 
                                edge_constraints=edge_constraints)
      #browser()
      return(list(value=res$value,init=res$init))
    }
    #browser()
    test0=grid_search(rho_vec[1])
    test1=grid_search(rho_vec[num_grid_sim])
    # parallelize process
    cores=detectCores()
    cl <- makeCluster(min(cores[1]-1,ncores)) #not to overload your computer
    registerDoParallel(cl)
    #browser()
    this <- foreach(i=1:length(rho_vec), .combine=cbind, 
                    .packages=PACKAGES, 
                    .export=c(names(as.list(.GlobalEnv)),ls())) %dopar% {
      grid_search(rho_vec[i])
    }
    stopCluster(cl)
    
    res = sapply((1:num_grid_sim), function(s) this[,s]$value)#numeric(num_grid_sim)
    
    par(mfrow = c(1, 1))
    plot(res,type="l",ylab="Translated Frob. norm",xlab="i")
    #browser()
    init <- c(this[,which.min(res)]$init,log(rho_vec[which.min(res)]))
    
    null_vec = which(round(c(1-sum(exp(init)[1:(length(exp(init))-1)]),
                             exp(init)[1:(length(exp(init))-1)]),15)<=0)
    if(length(null_vec)>0){
      matList$Gl[[1]] = calc_tilde_G_inv(matList$Ml[[1]],matList$Al[[1]],
                                         exp(init)[length(init)])[id_min,id_min]
      matList_full = link(matList)
      matList_full_extended = c(list(matrix(0,n,n)),matList_full)
      for(mat in matList_full_extended){
        diag(mat)=0
      }
      
      #null_vec = unique(c(1,null_vec))
      
      one_vec = (1:length(init))[-unique(c(1,null_vec))]#the null matrix can never be a target,
      # since its correlation is always 0
      
      ## choose vector pair with smallest distance in supports ##
      dist_matrix = sapply(null_vec, 
                           function(s) sapply(one_vec,
                                function(t) mat_support_distance(matList_full_extended[[s]],
                                                                 matList_full_extended[[t]])))
      
      #browser()
      if(is.vector(dist_matrix)){
        #if there is only one option, choose the one
        min_vec = sapply(seq_along(null_vec), function(s) which.min(dist_matrix[s]))
        arg_min1 = which.min(sapply(seq_along(null_vec), function(s) dist_matrix[s]))
        arg_min2 = 1
      } else{
        min_vec = sapply(seq_along(null_vec), function(s) which.min(dist_matrix[,s]))
        arg_min1 = which.min(sapply(seq_along(null_vec), function(s) dist_matrix[min_vec[s],s]))
        arg_min2 = min_vec[arg_min1]
      }
      #dist_matrix[arg_min2,arg_min1]
      r_min = null_vec[arg_min1]-1 # b chosen for the constraint of the form b>a/K
      s_min = one_vec[arg_min2]-1 # a chosen for the constraint of the form b>a/K
      #browser()
      constraint_digit = ((exp(init)[-length(init)])[s_min])*
        mean(matList_full_extended[[s_min+1]][matList_full_extended[[s_min+1]]>0])/(length(matList_full)+1)
      #browser()
      
      #In the case that ALL of the values are 0
      if(length(one_vec)==0){
        s_min=0
        r_min=1
        constraint_digit = 1e-15
      }
      
      edge_constraints[[counter]] = list(r_min=r_min, s_min=s_min, 
                                         constraint_digit=max(constraint_digit,1e-15))
      ## End: choose vector pair with smallest distance in supports ##
    } else{
      is_on_edge = FALSE
    }
  }
  #browser()
  return(init)
}

fit_cov = function(id_min, p, matList, init, Y=NULL,
                   link=function(matList) c(matList$Fk,matList$Gl),
                   link_der_rho=link_der_simple){
  #browser()
  # Optim
  n = dim(matList$Fk[[1]])[1]
  LogLikLogParm = function(x) LogLikLogParm_02(id_min=id_min, logParm=x, matList=matList, Y=Y, link=link)
  GradLogLikLogParm = function(x) GradLogLikLogParm_02(id_min=id_min, logParm=x, matList=matList, Y=Y, link=link, link_der_rho=link_der_rho)

  #browser()
  logLikInit <- LogLikLogParm(init)
  #browser()
  
  # fit3 <- try(optim(init, fn=LogLikLogParm, control=list(fnscale=-1, trace=1), method='Nelder-Mead'))
  # if(!is.character(fit3[1])){
  #   SigmaHat3 <- CovMat_03(id_min=id_min,parm=backward_transform_param(fit3$par), matList=matList, link=link)$Sigma
  #   param_fit3 = backward_transform_param(fit3$par)
  # } else{
  #   SigmaHat3 = NULL
  #   param_fit3 = NULL
  # }
  # #browser()
  # fit3 <- try(optim(init, fn=LogLikLogParm, control=list(fnscale=-1, trace=1), method='BFGS'))
  # if(!is.character(fit3[1])){
  #   SigmaHat3 <- CovMat_03(id_min=id_min, parm=backward_transform_param(fit3$par), matList=matList, link=link)$Sigma
  #   param_fit3 = backward_transform_param(fit3$par)
  # } else{
  #   SigmaHat3 = NULL
  #   param_fit3 = NULL
  # }
  #browser()
  fit3 <- try(optim(par=init, fn=LogLikLogParm, gr=GradLogLikLogParm, control=list(fnscale=-1, trace=1,maxit=500), method='BFGS'))
  if(!is.character(fit3[1])){
    SigmaHat3 <- CovMat_03(id_min=id_min,parm=backward_transform_param(fit3$par), matList=matList, link=link)$Sigma
    param_fit3 = backward_transform_param(fit3$par)
  } else{
    SigmaHat3 = NULL
    param_fit3 = NULL
  }
  #browser()
  bic = -2*true_LogLikParm_02(id_min, param_fit3, matList, Y, link=link) + length(init)*log(dim(Y)[1])
  #bic = -2*fit3$value + length(init)*log(dim(Y)[1])
  #browser()
  return(list(SigmaHat1=SigmaHat3,#SigmaHat2=SigmaHat2,SigmaHat3=SigmaHat3,
              param_fit1=param_fit3,# param_fit2=param_fit2, param_fit3=param_fit3,
              bic=bic))
}
