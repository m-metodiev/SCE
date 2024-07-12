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