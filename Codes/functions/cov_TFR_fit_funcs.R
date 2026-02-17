PACKAGES_MISC = c("sna", "xts") # stuff
PACKAGES_STAT_OPTIM = c("mvtnorm","quadprog","pracma", "Matrix","missMDA", "covFactorModel",
                        "quadprog", "cvCovEst", "CovTools", "CVglasso") # stat,math stuff
#PACKAGES_VIS = c("UpSetR", "ggplot2", "grid", "plyr", "igraph") # plot stuff
PACKAGES_PARALLEL = c("foreach", 
                      "parallel", "doParallel") # parallelization stuff
PACKAGES = c(PACKAGES_MISC, PACKAGES_STAT_OPTIM, 
             PACKAGES_PARALLEL) # all of the stuff
lapply(PACKAGES, require, character.only = TRUE)

fit_param = function(dimension, num_observations, matList, sim, id_min,
                     filename_param_fit=NULL, filename_ests=NULL, filename_bic=NULL, 
                     filename_error_measures=NULL,
                     num_grid_sim = 100, ncores=8,
                     save=TRUE,Sigma=NULL,
                     link=function(matList) c(matList$Fk,matList$Gl), 
                     link_der_beta=link_der_simple,
                     normalize_data=FALSE, compute_WSCE=FALSE,
                     simple=FALSE, init=NULL){
  
  matList$Gl[[1]] = matrix(0,dimension,dimension) # Just a dummy-variable, will be declared later
  num_observations_is_one=FALSE
  if(num_observations==1){
    sim$Y = t(matrix(sim$Y))
    num_observations_is_one=TRUE
  } # code is slightly different if num_observations==1

  # the SCE needs to be computed on the normalized dataset
  if(normalize_data){
    simY_normalized = (sim$Y - t(matrix(rep(colMeans(sim$Y),num_observations),ncol=num_observations)))%*%diag(1/apply(sim$Y,2,sd))
    sim$corY=cor(sim$Y)
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

  exp_param_optim_frob = exp(init)
  Sigma_0_optim_frob = CovMat_03(parm=exp_param_optim_frob, matList,id_min=id_min)$Sigma

  list2env(fit_cov(id_min, num_observations,
                      matList,
                      forward_transform_param(exp(init)),
                      as.matrix(simY_normalized),
                        link=link,link_der_beta=link_der_beta), envir = .GlobalEnv)

  if(sum(is.na(sim$Y))>0){
    has_missingvalues=TRUE
    Y_nonmissing = imputePCA(sim$Y,ncp=length(init))$completeObs
  } else{
    has_missingvalues=FALSE
    Y_nonmissing = sim$Y
  }
  
  if((num_observations>1) & (!simple)){
    ests = list(as.matrix(cor(Y_nonmissing)),
                cov2cor(covFactorModel::covFactorModel(xts(x=Y_nonmissing, order.by=Sys.Date()-1:num_observations),K=length(init))),
                cov2cor(CVglasso(X=Y_nonmissing)$Sigma),
                cov2cor(linearShrinkLWEst(Y_nonmissing)),
                as.matrix(Sigma_0_optim_frob),
                as.matrix(SigmaHat1))
    
    if(compute_WSCE){
      if(normalize_data){
        pearson_mat = cor(sim$Y) # mean and variance unknown
        use_bootstrap = TRUE
      } else{
        pearson_mat = cor_from_standard_errors(Y_nonmissing) # mean and variance known
        use_bootstrap = FALSE # used anyway if part of the data is missing
      }
      lambda = compute_lambda_opt(id_min,parm=param_fit1,matList=matList,
                                  link=link,link_der_beta=link_der_beta,Y=sim$Y,
                                  pearson_mat=pearson_mat,SCE_mat=as.matrix(SigmaHat1),
                                  Y_nonmissing=Y_nonmissing, use_bootstrap=use_bootstrap)
      print("lambda")
      print(lambda)
      ests = c(ests, list(lambda*as.matrix(SigmaHat1)+(1-lambda)*as.matrix(cor(Y_nonmissing))))
    }
  } else{
    ests = list(as.matrix(Sigma_0_optim_frob),as.matrix(SigmaHat1))
  } # If num_observations=1, most estimators cant be computed  

  #browser()
  
  if(save==TRUE){
    if(length(param_fit1)==7){
      param_fit1=c(param_fit1[1:4],0,param_fit1[5:7])
    }
    write_param(t(c(param_fit1)),filename=filename_param_fit,
                matList = matList, link=link)
    write.csv(sapply(seq_along(ests),function(s) ests[[s]]),file=filename_ests)
    write.csv(as.data.frame(bic),file=filename_bic)
    if(!is.null(filename_error_measures)){
      write_summary_measures(filename_ests = filename_ests, 
                             filename_error_measures = filename_error_measures,
                             num_observations_is_one=num_observations_is_one,Sigma=Sigma)
    }
  } else{
    if(!is.null(filename_error_measures)){
      summary_measures = sapply(1:length(ests),
                                function(i) c(mean(abs(ests[[i]]-Sigma)),
                                              sqrt(mean((ests[[i]]-Sigma)^2))))
      if(compute_WSCE){
        return(list(summary_measures, bic, ests, lambda, param_fit1))
      } else{
        return(list(summary_measures, bic, ests, param_fit1))
      }
    } else{ # error measures aren't necessarily always available
      return(list(t(c(param_fit1)), ests, bic))
    }
  }
  #browser()
  
}

fit_init = function(id_min, matList, corY,
                    num_grid_sim = 100, ncores=8,
                    link=function(matList) c(matList$Fk,matList$Gl)){

  dimension = dim(matList$Fk[[1]])[1]
  s = dim(matList$Al[[1]])[1]
  #Need grid-search for beta because the norm is not quadratic w.r.t. beta
  beta=1.5
  xi = (1:(num_grid_sim+1))/(num_grid_sim+1)
  #tan-hyperbolic-spaced grid because beta approaches 1
  beta_vec = (1-tanh(beta*(1+xi))/tanh(beta))/(min((1-tanh(beta*(1+xi))/tanh(beta))))
  beta_vec = beta_vec[-length(beta_vec)]
  
  is_on_edge = TRUE # solution can lie on the edge of the parameter space
  edge_constraints = list() # params which lie on the edge will be adjusted in constraints
  
  counter = 0
  while(is_on_edge){
    counter = counter + 1
    grid_search = function(beta){
      matList$Gl[[1]] = calc_tilde_G_inv(M=matList$Ml[[1]],A=matList$Al[[1]],beta=beta,
                                         U_full=matList$U_full, solve_U_full=matList$solve_U_full,
                                         solve_M_no_islands=matList$solve_M_no_islands,
                                         eigen_real=matList$eigen_real)[id_min,id_min]
      matList_full = link(matList)
      res = calc_Sigma_opt_frob(matList_full, corY, 
                                edge_constraints=edge_constraints)
      return(list(value=res$value,init=res$init))
    }

    test0=grid_search(beta_vec[1])
    test1=grid_search(beta_vec[num_grid_sim])
    # parallelize process
    cores=detectCores()
    cl <- makeCluster(min(cores[1]-1,ncores)) # to not overload your computer
    registerDoParallel(cl)
    this <- foreach(i=1:length(beta_vec), .combine=cbind, 
                    .packages=PACKAGES, 
                    .export=c(names(as.list(.GlobalEnv)),ls())) %dopar% {
      grid_search(beta_vec[i])
    }
    stopCluster(cl)
    
    res = sapply((1:num_grid_sim), function(s) this[,s]$value)
    
    par(mfrow = c(1, 1))
    plot(res,type="l",ylab="Translated Frob. norm",xlab="i")

    init <- c(this[,which.min(res)]$init,log(beta_vec[which.min(res)]))
    
    null_vec = which(round(c(1-sum(exp(init)[1:(length(exp(init))-1)]),
                             exp(init)[1:(length(exp(init))-1)]),15)<=0)
    if(length(null_vec)>0){
      matList$Gl[[1]] = calc_tilde_G_inv(matList$Ml[[1]],matList$Al[[1]],
                                         exp(init)[length(init)])[id_min,id_min]
      matList_full = link(matList)
      matList_full_extended = c(list(matrix(0,dimension,dimension)),matList_full)
      for(mat in matList_full_extended){
        diag(mat)=0
      }
      
      one_vec = (1:length(init))[-unique(c(1,null_vec))]#the null matrix can never be a target,
      # since its correlation is always 0
      
      ## choose vector pair with smallest distance in supports ##
      dist_matrix = sapply(null_vec, 
                           function(s) sapply(one_vec,
                                function(t) mat_support_distance(matList_full_extended[[s]],
                                                                 matList_full_extended[[t]])))

      if(is.vector(dist_matrix)){
        # if there is only one option, choose the one
        min_vec = sapply(seq_along(null_vec), function(s) which.min(dist_matrix[s]))
        arg_min1 = which.min(sapply(seq_along(null_vec), function(s) dist_matrix[s]))
        arg_min2 = 1
      } else{
        min_vec = sapply(seq_along(null_vec), function(s) which.min(dist_matrix[,s]))
        arg_min1 = which.min(sapply(seq_along(null_vec), function(s) dist_matrix[min_vec[s],s]))
        arg_min2 = min_vec[arg_min1]
      }
      r_min = null_vec[arg_min1]-1 # b chosen for the constraint of the form b>a/K
      s_min = one_vec[arg_min2]-1 # a chosen for the constraint of the form b>a/K
      constraint_digit = ((exp(init)[-length(init)])[s_min])*
        mean(matList_full_extended[[s_min+1]][matList_full_extended[[s_min+1]]>0])/(length(matList_full)+1)
      # K chosen for the constraint of the form b>a/K
      
      # In the case that ALL of the values are 0
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
  return(init)
}

#' Title
#'
#' @param id_min 
#' @param num_observations 
#' @param matList 
#' @param init 
#' @param Y 
#' @param link 
#' @param link_der_beta 
#'
#' @return
#' @export
#'
#' @examples
fit_cov = function(id_min, num_observations, matList, init, Y=NULL,
                   link=function(matList) c(matList$Fk,matList$Gl),
                   link_der_beta=link_der_simple){

  n = dim(matList$Fk[[1]])[1]
  LogLikLogParm = function(x) LogLikLogParm_02(id_min=id_min, logParm=x, matList=matList, Y=Y, link=link)
  GradLogLikLogParm = function(x) GradLogLikLogParm_02(id_min=id_min, logParm=x, matList=matList, Y=Y, link=link, link_der_beta=link_der_beta)

  logLikInit <- LogLikLogParm(init)

  fit3 <- try(optim(par=init, fn=LogLikLogParm, gr=GradLogLikLogParm, control=list(fnscale=-1, trace=1,maxit=500), method='BFGS'))
  if(!is.character(fit3[1])){
    SigmaHat3 <- CovMat_03(id_min=id_min,parm=backward_transform_param(fit3$par), matList=matList, link=link)$Sigma
    param_fit3 = backward_transform_param(fit3$par)
  } else{
    SigmaHat3 = NULL
    param_fit3 = NULL
  }
  # fit3 <- try(optim(init, fn=LogLikLogParm, control=list(fnscale=-1, trace=1), method='BFGS'))
  # if(!is.character(fit3[1])){
  #   SigmaHat3 <- CovMat_03(id_min=id_min, parm=backward_transform_param(fit3$par), matList=matList, link=link)$Sigma
  #   param_fit3 = backward_transform_param(fit3$par)
  # } else{
  #   SigmaHat3 = NULL
  #   param_fit3 = NULL
  # }
  print("INIT:")
  print(backward_transform_param(init))
  bic = -2*true_LogLikParm_02(id_min, param_fit3, matList, Y, link=link) + length(init)*log(dim(Y)[1])

  return(list(SigmaHat1=SigmaHat3,
              param_fit1=param_fit3,
              bic=bic))
}
