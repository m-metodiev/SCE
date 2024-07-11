rm(list=ls())
source("functions/FunctionsCovTFR_02.R")
source("functions/funcs_frob_norm_opt.R")
source("functions/cov_TFR_sim_funcs.R")
source("functions/cov_TFR_fit_funcs.R")
source("functions/cov_TFR_plot_funcs.R")
source("functions/cov_TFR_data_funcs.R")

seed <- 3; set.seed(seed)
data_source = "data/"

#### Simulation 1: Simulate both, the matrices and the standardized error ####

### Initialization ###

n <- 200; p <- 10; rho = .982

# Parms
#alpha <- c(.11,.05,.09)#beta <- c(.01)#

matList2 = sim_matList(n,rho=rho,num_F=2,k_vec=c(3,10),num_G=1,F_0=FALSE)
matList2$Fk[[length(matList2$Fk)+1]] = matrix(1,nrow=n,ncol=n)

write_matList(matList=matList2, filename=paste(data_source,"sim_01_matList.csv",sep=""))


parm = c(.05,0.09,.11,.74,rho)

#true_param_mean_neighbor_effect = calc_mean_neighbor_effect(matList2,parm[5],parm[4],1:dim(matList2$Fk[[1]])[1])
write_param(t(parm),matList=matList2,
            filename=paste(data_source,"sim_01_true_param.csv",sep=""))

#check that the maximum neighbor effect is around .26
test = calc_tilde_G_inv(matList2$Ml[[1]],matList2$Al[[1]],rho)
diag(test)=0
max(test)*.74#maximum neighbor effect around .26
mean(test[matList2$Al[[1]]!=0])*.74#neighbor effect around .21

Sigma = CovMat_03(parm, matList2)$Sigma

sim2 = sim_cov(p, Sigma)

source("functions/cov_TFR_fit_funcs.R")
fit_param(n, p, matList2, sim2, Sigma=Sigma,id_min=1:dim(matList2$Fk[[1]])[1],
          filename_param_fit=paste(data_source,"sim_01_param_fit.csv",sep=""),
          filename_ests = paste(data_source,"sim_01_ests.csv",sep=""),
          filename_bic = paste(data_source,"sim_01_bic.csv",sep=""),
          filename_error_measures = paste(data_source,"sim_01_error_measures.csv",sep=""),
          compute_WSCE = TRUE)

read.csv(file=paste(data_source,"sim_01_error_measures.csv",sep=""))

fit_param(n, p, matList2, sim2, Sigma=Sigma,id_min=1:dim(matList2$Fk[[1]])[1],
          filename_param_fit=paste(data_source,"sim_01_param_fit_musigma_unknown.csv",sep=""),
          filename_ests = paste(data_source,"sim_01_ests_musigma_unknown.csv",sep=""),
          filename_bic = paste(data_source,"sim_01_bic_musigma_unknown.csv",sep=""),
          filename_error_measures = paste(data_source,"sim_01_error_measures_musigma_unknown.csv",sep=""),
          normalize_data = TRUE, compute_WSCE = TRUE)
read.csv(file=paste(data_source,"sim_01_error_measures_musigma_unknown.csv",sep=""))

## Compute performance over repeated simulations
set.seed(seed)
num_sim=10
sim_func = function(seed) sim_cov_with_seed(p, Sigma, seed)
fit1_func = function(sim) fit_param(n, p, matList2, sim, id_min=1:dim(matList2$Fk[[1]])[1],
                                    filename_error_measures=TRUE,save=FALSE,Sigma=Sigma, compute_WSCE = TRUE)
fit2_func = function(sim) fit_param(n, p, matList2, sim, id_min=1:dim(matList2$Fk[[1]])[1],
                                    filename_error_measures = TRUE,
                                    save=FALSE,Sigma=Sigma,link=combined_matList, link_der_rho = link_der_combined,
                                    compute_WSCE=TRUE)
sims_errors_and_bic = sim_errors_and_bic(sim_func,fit1_func,fit2_func,num_sim)
write.csv(sims_errors_and_bic,
          file=paste(data_source,"sim_01_sims_errors_and_bic.csv",sep=""))

set.seed(seed)
num_sim=10
sim_func = function(seed) sim_cov_with_seed(p, Sigma, seed)
fit1_func = function(sim) fit_param(n, p, matList2, sim, id_min=1:dim(matList2$Fk[[1]])[1],
                                    filename_error_measures=TRUE,save=FALSE,Sigma=Sigma,
                                    compute_WSCE = TRUE, normalize_data = TRUE)
fit2_func = function(sim) fit_param(n, p, matList2, sim, id_min=1:dim(matList2$Fk[[1]])[1],
                                    filename_error_measures = TRUE,
                                    save=FALSE,Sigma=Sigma,link=combined_matList, link_der_rho = link_der_combined,
                                    compute_WSCE = TRUE, normalize_data = TRUE)
sims_errors_and_bic = sim_errors_and_bic(sim_func,fit1_func,fit2_func,num_sim)
write.csv(sims_errors_and_bic,
          file=paste(data_source,"sim_01_sims_errors_and_bic_musigma_unknown.csv",sep=""))

# The following code is not used anymore, since this simulation is now done
# in the setting in which the F_k and G matrices are non-random.
# set.seed(seed)
# n_vec = c(10, 50, 100, 200, 500, 1000)
# res_list_list = list()
# 
# #The "full" sample with sample size n=1000, but T=1
# set.seed(seed)
# n = 1000; p <- 1; rho = .982
# matList0 = sim_matList(n,rho=rho,num_F=2,k_vec=c(3,10),num_G=1,F_0=FALSE)
# matList0$Fk[[length(matList0$Fk)+1]] = matrix(1,nrow=n,ncol=n)
# parm = c(.05,0.09,.11,.74,rho)
# Sigma = CovMat_03(parm, matList0)$Sigma
# 
# car_matrices = calc_tilde_G_inv(M=matList0$Ml[[1]],A=matList0$Al[[1]],rho=rho, return_U_D_M = TRUE)
# matList0 = c(matList0,car_matrices)
# 
# for(N in 1:10){
#   res_list = list()
#   sim0 = sim_cov(p, Sigma)
#   sim0$corY = (((t(sim0$Y)%*%sim0$Y)>0)-((t(sim0$Y)%*%sim0$Y)<0)) + 0
#   sim0$corY = t(sim0$Y)%*%sim0$Y
#   
#   matList1 = matList0
#   sim1 = sim0
#   
#   source("cov_TFR_fit_funcs.R")
#   for(i in seq_along(n_vec)){
#     #browser()
#     n <- n_vec[i];# p <- 1; rho = .982
#     for(j in seq_along(matList1$Fk)){
#       matList1$Fk[[j]] = matList0$Fk[[j]][1:n,1:n]
#     } # truncate the Fk, because we only need this part of the info
#     
#     sim1$Y = sim0$Y[1,1:n]
#     sim1$corY = sim0$corY[1:n,1:n]
#     Sigma_1 = Sigma[1:n,1:n]
#     
#     res_list[[i]] = fit_param(n, p, matList1, sim1,
#                               id_min=1:dim(matList1$Fk[[1]])[1], save=FALSE,
#                               filename_error_measures=TRUE,Sigma=Sigma_1)
#     #browser()
#   }
#   res_list_list[[N]] = res_list
# }
# 
# save(res_list_list,file="test0002.Rdata")
# length(res_list_list)
# res_matrix = sapply(1:length(res_list_list), function (t) 
#   sapply(1:length(n_vec), function(s) res_list_list[[t]][[s]][[1]][2,2]))
# write.csv(res_matrix,file=paste(data_source,"sim_00_sample_one.csv",sep=""))