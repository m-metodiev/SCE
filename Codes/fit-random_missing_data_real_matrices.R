rm(list=ls())
source("functions/FunctionsCovTFR_02.R")
source("functions/funcs_frob_norm_opt.R")
source("functions/cov_TFR_sim_funcs.R")
source("functions/cov_TFR_fit_funcs.R")
source("functions/cov_TFR_plot_funcs.R")
source("functions/cov_TFR_data_funcs.R")

seed <- 3; set.seed(seed)
data_source = "data/"

#### Simulation 2: Matrices from the real data ####

#### TODO ####

#Create new mean
# seed <- 3; set.seed(seed)
# set.seed(seed)
# n <- 195; p <- 1; rho = .982
# matList0 = sim_matList(n,rho=rho,num_F=2,k_vec=c(3,10),num_G=1,F_0=FALSE)
# matList0$Fk[[length(matList0$Fk)+1]] = matrix(1,nrow=n,ncol=n)
# parm = c(.05,0.09,.11,.74,rho)
# Sigma2 = CovMat_03(parm, matList0)$Sigma

#### TODO ####

seed=3
set.seed(seed)

source("functions/get_data_covar.R")

read_plot = read_plot_FITcomps_std(filename="../Data/TFR_pieces_202311/standardized_residuals_202311/FITcomps_std_residuals_sample%i_202311.txt") # Initializing and plotting the real data 
FITcomps_std_total = read_plot$FITcomps_std_total
FITcomps_std = read_plot$FITcomps_std

read_names = read_names_FITcomps_std_total(FITcomps_std_total, covar, model="all_values")
names_by_id = read_names$names_by_id
all_min = read_names$all_min

preproc_res = preproc_FITcomps_std(all_min, names_by_id, FITcomps_std, covar)
matList_final = preproc_res$matList_final
dim(matList_final$Fk[[1]])
id_min = preproc_res$id_min

# Parms
n <- 195; p <- 11; rho = .35
#alpha <- c(.11,.05,.09)#beta <- c(.01)#

matList2 = matList_final#sim_matList(n,rho=rho,num_F=2,k_vec=c(3,10),num_G=1,F_0=FALSE)
write_matList(matList=matList2, filename=paste(data_source,"sim_02_matList.csv",sep=""))


#matList2$Fk[[length(matList2$Fk)+1]] = matrix(1,nrow=n,ncol=n)

parm = c(.05,0.09,.11,.74,rho)
write_param(t(parm),matList=matList2,filename=paste(data_source,"sim_03_true_param.csv",sep=""))

test = calc_tilde_G_inv(matList2$Ml[[1]],matList2$Al[[1]],rho)[id_min,id_min]
diag(test)=0
max(test)*.74#maximum neighbor effect around .26
G_inv = calc_tilde_G_inv(matList2$Ml[[1]],matList2$Al[[1]],rho)[id_min,id_min]
A = matList2$Al[[1]][id_min,id_min]
A[is.na(A)] = 0 # since islands can return NA-values

covMat <- CovMat_03(parm, matList2,id_min=id_min)
Sigma <- covMat$Sigma#cov2cor(covMat$Sigma)
sim2 = sim_cov(p, as.matrix(Sigma))

#### Simulation 3: Matrices and missing values from the real data ####
Y = FITcomps_std_total[2:12,which(all_min==1)]
Y = Y[,sapply(id_min,function(id) which(preproc_res$FITcomps_std_iso==preproc_res$iso_id_key[id]))]

sim3 = sim_cov(p, as.matrix(Sigma))
sim3$Y[is.na(Y)] = NA

corY = matrix(ncol=n,nrow=n)

#use pairwise correlation estimates
for(i in (1:n)){
  for(j in (i:n)){
    #browser()
    Yi_notmissing = !is.na(sim3$Y[,i])
    Yj_notmissing = !is.na(sim3$Y[,j])
    corY_ij = cor(sim3$Y[Yi_notmissing&Yj_notmissing,i], sim3$Y[Yi_notmissing&Yj_notmissing,j])
    corY[i,j] = corY_ij
    corY[j,i] = corY_ij
  }
}


sim3$corY = corY
matList3 = matList2

source("functions/cov_TFR_fit_funcs.R")
fit_param(n, p, matList3, sim3, id_min=id_min, Sigma=Sigma,
          filename_param_fit=paste(data_source,"sim_03_param_fit.csv",sep=""),
          filename_ests = paste(data_source,"sim_03_ests.csv",sep=""),
          filename_bic = paste(data_source,"sim_03_bic.csv",sep=""),
          filename_error_measures = paste(data_source,"sim_03_error_measures.csv",sep=""),
          compute_WSCE = TRUE)

read.csv(file=paste(data_source,"sim_03_error_measures.csv",sep=""))

fit_param(n, p, matList3, sim3, id_min=id_min, Sigma=Sigma,
          filename_param_fit=paste(data_source,"sim_03_param_fit_combined_effects.csv",sep=""),
          filename_ests = paste(data_source,"sim_03_ests_combined_effects.csv",sep=""), 
          filename_bic = paste(data_source,"sim_03_bic_combined_effects.csv",sep=""),
          filename_error_measures = paste(data_source,"sim_03_error_measures_combined_effects.csv",sep=""),
          link=combined_matList)
read.csv(file=paste(data_source,"sim_03_error_measures_combined_effects.csv",sep=""))

read.csv(file=paste(data_source,"sim_03_bic_combined_effects.csv",sep=""))$bic - 
  read.csv(file=paste(data_source,"sim_03_bic.csv",sep=""))$bic

set.seed(seed)

num_sim=10
sim_func = function(seed){
  set.seed(seed)
  sim_covs = sim_cov(p, as.matrix(Sigma))
  sim_covs$Y[is.na(Y)] = NA
  corY = matrix(ncol=n,nrow=n)
  
  #use pairwise correlation estimates
  for(i in (1:n)){
    for(j in (i:n)){
      #browser()
      Yi_notmissing = !is.na(sim_covs$Y[,i])
      Yj_notmissing = !is.na(sim_covs$Y[,j])
      corY_ij = cor(sim_covs$Y[Yi_notmissing&Yj_notmissing,i], sim_covs$Y[Yi_notmissing&Yj_notmissing,j])
      corY[i,j] = corY_ij
      corY[j,i] = corY_ij
    }
  }
  sim_covs$corY = corY
  return(sim_covs)
} 
fit1_func = function(sim) fit_param(n, p, matList2, sim, id_min=id_min,
                                    filename_error_measures=TRUE,save=FALSE,
                                    Sigma=Sigma, compute_WSCE = TRUE)
fit2_func = function(sim) fit_param(n, p, matList2, sim, id_min=id_min,
                                    filename_error_measures = TRUE,
                                    link=combined_matList, save=FALSE,
                                    Sigma=Sigma)
source("functions/cov_TFR_fit_funcs.R")
sims_errors_and_bic = sim_errors_and_bic(sim_func,fit1_func,fit2_func,num_sim)
write.csv(sims_errors_and_bic, file=paste(data_source,"sim_03_sims_errors_and_bic.csv",sep=""))