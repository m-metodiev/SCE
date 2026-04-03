PACKAGES_MISC = c("sna", "xts") # stuff
PACKAGES_STAT_OPTIM = c("mvtnorm","quadprog","pracma", "Matrix","missMDA", "covFactorModel",
                        "quadprog", "cvCovEst", "CovTools", "CVglasso") # stat,math stuff
#PACKAGES_VIS = c("UpSetR", "ggplot2", "grid", "plyr", "igraph") # plot stuff
PACKAGES_PARALLEL = c("foreach", 
                      "parallel", "doParallel") # parallelization stuff
PACKAGES = c(PACKAGES_MISC, PACKAGES_STAT_OPTIM, 
             PACKAGES_PARALLEL) # all of the stuff
lapply(PACKAGES, require, character.only = TRUE)
library(scov)

# if not positive semidefinite, the matrix is mapped to the
# positive definite matrix closest to it by using nearPD
fit_param = function(dimension, num_observations, matList, sim, adj_positions,
                     filename_param_fit=NULL, filename_ests=NULL, filename_bic=NULL, 
                     filename_error_measures=NULL,
                     grid_size = 100, ncores=5,
                     save=TRUE,Sigma=NULL,
                     #link=function(matList) list(matList_full=c(matList$Fk,matList$Gl)), 
                     #link_der_beta=link_der_simple,
                     normalize_data=FALSE, compute_WSCE=FALSE,
                     simple=FALSE, init=NULL, adj_matrix=NULL,
                     interaction_effects=list(),
                     num_bootstrap_iters=100,
                     parallelize = TRUE,
                     joint_estimation = FALSE,
                     already_demeaned = FALSE){

  if(is.null(adj_matrix)&(!is.null(matList$vois))){
    Ml_Al = get_M_A(adj_matrix=matList$vois)
    matList$Ml = Ml_Al$Ml
    matList$Al = Ml_Al$Al
    matList$Gl = NULL
  }
  
  matList$Gl[[1]] = matrix(0,dimension,dimension) # Just a dummy-variable, will be declared later
  num_observations_is_one=FALSE
  if(num_observations==1){
    sim$dataset = t(matrix(sim$dataset))
    num_observations_is_one=TRUE
  } # code is slightly different if num_observations==1

  # the SCE needs to be computed on the normalized dataset
  if(sum(is.na(sim$dataset))==0){
    if(normalize_data){
      mean_estim = colMeans(sim$dataset)
      sd_estim = apply(sim$dataset,2,sd)
      sim$correlation_matrix=cor(sim$dataset)
    } else{
      #simY_normalized = sim$Y
      if(!already_demeaned){
        mean_estim = rep(2.1,dimension)
      } else{
        mean_estim = rep(0,dimension)
      }
      
      sd_estim = rep(1,dimension)
      sim$correlation_matrix = cor_from_standard_errors(sim$dataset)
    }
    dataset_normalized = (sim$dataset - t(matrix(rep(mean_estim,num_observations),ncol=num_observations)))%*%diag(1/sd_estim)
  } else{
    # we do not (so far) simulate the case where the data is missing and
    # the mean and sd are unknown, so we assume that they are known
    if(!already_demeaned){
      mean_estim = rep(2.1,dimension)
    } else{
      mean_estim = rep(0,dimension)
    }
    sd_estim = rep(1,dimension)
    dataset_normalized = sim$dataset
  }
  # for the real dataset, we know that the standardized errors are 
  # already normalized, so we don't need to do that again

  if(is.null(init)){
    ive_estim = scov(pairwise_covariate_matrices = matList$Fk,
                     adj_matrix = matList$vois,
                     dataset = sim$dataset,
                     mean_estim = mean_estim, sd_estim = sd_estim,
                     grid_size = grid_size, ncores = ncores,
                     adj_positions=adj_positions,
                     interaction_effects=interaction_effects,
                     parallelize = parallelize,
                     semiparametric = TRUE)
    init = ive_estim$parm
  }
  
  Sigma_0_optim_frob = ive_estim$corrmat_estim
  #Sigma_0_optim_frob = CovMat_03(parm=init, matList,
  #                               adj_positions=adj_positions,
  #                               interaction_effects = interaction_effects)$Sigma
  sce_estim = scov(pairwise_covariate_matrices=matList$Fk,
                   adj_matrix = matList$vois,
                   dataset = sim$dataset,
                   mean_estim = mean_estim, sd_estim = sd_estim,
                   grid_size = grid_size, ncores = ncores,
                   adj_positions = adj_positions,
                   init = init,
                   interaction_effects=interaction_effects,
                   parallelize = parallelize,
                   semiparametric = FALSE,
                   misspecification = FALSE,
                   joint_estimation=joint_estimation)
  
  SigmaHat1 = sce_estim$corrmat_estim
  param_fit1 = sce_estim$parm
  bic = sce_estim$bic
  avg_effects = sce_estim$average_effects
  
  if(sum(is.na(sim$dataset))>0){
    has_missingvalues=TRUE
    dataset_nonmissing = imputePCA(sim$dataset,ncp=length(init))$completeObs
  } else{
    has_missingvalues=FALSE
    dataset_nonmissing = sim$dataset
  }
  
  if((num_observations>1) & (!simple)){
    ests = list(as.matrix(cor(dataset_nonmissing)),
                cov2cor(covFactorModel::covFactorModel(xts(x=dataset_nonmissing, order.by=Sys.Date()-1:num_observations),K=length(init))),
                cov2cor(CVglasso(X=dataset_nonmissing)$Sigma),
                cov2cor(linearShrinkLWEst(dataset_nonmissing)),
                as.matrix(Sigma_0_optim_frob),
                as.matrix(SigmaHat1))
    
    if(compute_WSCE){
      if(normalize_data){
        pearson_mat = cor(sim$dataset) # mean and variance unknown
        use_bootstrap = TRUE
      } else{
        if(!already_demeaned){
          pearson_mat = cor_from_standard_errors(dataset_nonmissing - 2.1)
        } else{
          pearson_mat = cor_from_standard_errors(dataset_nonmissing)
        }
        # mean and variance known
        use_bootstrap = FALSE # used anyway if part of the data is missing
      }
      wsce_estim = scov(pairwise_covariate_matrices=matList$Fk,
                        adj_matrix=matList$vois,
                        dataset=sim$dataset,
                        mean_estim = mean_estim, sd_estim = sd_estim,
                        grid_size=grid_size, parallelize = parallelize, 
                        ncores=ncores, adj_positions = adj_positions,
                        interaction_effects = interaction_effects,
                        init = init,
                        use_bootstrap=use_bootstrap,
                        num_bootstrap_iters=num_bootstrap_iters,
                        semiparametric=FALSE,
                        misspecification=TRUE,
                        joint_estimation=joint_estimation)
      
      lambda = wsce_estim$lambda
      print("lambda")
      print(lambda)
      ests = c(ests, list(wsce_estim$corrmat_estim))
    }
  } else{
    ests = list(as.matrix(Sigma_0_optim_frob),as.matrix(SigmaHat1))
  } # If num_observations=1, most estimators cant be computed  
  
  if(save==TRUE){
    if(length(param_fit1)==7){
      param_fit1=c(param_fit1[1:4],0,param_fit1[5:7])
    }
    write_param(t(c(param_fit1)),filename=filename_param_fit,
                matList = matList, interaction_effects=interaction_effects)
    write.csv(sapply(seq_along(ests),function(s) ests[[s]]),file=filename_ests)
    write.csv(as.data.frame(bic),file=filename_bic)
    if(!is.null(filename_error_measures)){
      write_summary_measures(filename_ests = filename_ests, 
                             filename_error_measures = filename_error_measures,
                             Sigma=Sigma)
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
      return(list(t(c(param_fit1)), ests, bic, avg_effects))
    }
  }
  # 
  
}

#' Computes non-positive-semidefinite approximation of correlation matrix
#'
#' @param dataset n x d matrix (n = number of observations, d = dimension)
#'
#' @returns non-positive-semidefinite approximation of correlation matrix
#' @keywords internal
compute_marginal_cor = function(dataset){
  dimension = ncol(dataset)
  correlation_matrix = matrix(ncol=dimension,nrow=dimension)
  
  #use pairwise correlation estimates
  for(i in (1:dimension)){
    for(j in (i:dimension)){
      dataseti_notmissing = !is.na(dataset[,i])
      datasetj_notmissing = !is.na(dataset[,j])
      correlation_matrix_ij = cor(dataset[dataseti_notmissing&
                                            datasetj_notmissing,i],
                                  dataset[dataseti_notmissing&
                                            datasetj_notmissing,j])
      correlation_matrix[i,j] = correlation_matrix_ij
      correlation_matrix[j,i] = correlation_matrix_ij
    }
  }
  return(correlation_matrix)
}