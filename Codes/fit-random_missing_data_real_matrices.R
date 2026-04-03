rm(list=ls())
source("functions/FunctionsCovTFR_02.R")
source("functions/funcs_frob_norm_opt.R")
source("functions/cov_TFR_sim_funcs.R")
source("functions/cov_TFR_fit_funcs.R")
source("functions/cov_TFR_plot_funcs.R")
source("functions/cov_TFR_data_funcs.R")

seed <- 3; set.seed(seed)
data_source = "data/"

#### Simulation 3: Matrices from the real data with missing values ####

## BEGIN Initialization ##
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
adj_positions = preproc_res$adj_positions

# Parms
dimension <- 195; num_observations <- 11; beta = .35

matList2 = matList_final
write_matList(matList=matList2, filename=paste(data_source,"sim_02_matList.csv",sep=""))

parm = c(.05,0.09,.11,.74,beta)
# write_param(t(parm),matList=matList2,filename=paste(data_source,"sim_03_true_param.csv",sep=""))

test = tilde_G_inv(matList2$Ml[[1]],matList2$Al[[1]],beta)[adj_positions,adj_positions]
diag(test)=0
max(test)*.74 # maximum neighbor effect around .26
G_inv = tilde_G_inv(matList2$Ml[[1]],matList2$Al[[1]],beta)[adj_positions,adj_positions]
A = matList2$Al[[1]][adj_positions,adj_positions]
A[is.na(A)] = 0 # since islands can return NA-values

matList2$Gl[[1]] = G_inv
covMat <- CovMat_03(parm, matList2,adj_positions=adj_positions)
Sigma <- covMat$Sigma
sim2 = sim_cov(num_observations, as.matrix(Sigma))

dataset = FITcomps_std_total[2:12,which(all_min==1)]
dataset = dataset[,sapply(adj_positions,function(id) which(preproc_res$FITcomps_std_iso==preproc_res$iso_id_key[id]))]

sim3 = sim_cov(num_observations, as.matrix(Sigma))
sim3$dataset[is.na(dataset)] = NA

correlation_matrix = compute_marginal_cor(dataset)
# correlation_matrix = matrix(ncol=dimension,nrow=dimension)
# 
# #use pairwise correlation estimates
# for(i in (1:dimension)){
#   for(j in (i:dimension)){
#     #browser()
#     dataseti_notmissing = !is.na(sim3$dataset[,i])
#     datasetj_notmissing = !is.na(sim3$dataset[,j])
#     correlation_matrix_ij = cor(sim3$dataset[dataseti_notmissing&datasetj_notmissing,i], 
#                   sim3$dataset[dataseti_notmissing&datasetj_notmissing,j])
#     correlation_matrix[i,j] = correlation_matrix_ij
#     correlation_matrix[j,i] = correlation_matrix_ij
#   }
# }


sim3$correlation_matrix = correlation_matrix
matList3 = matList2
names(matList3$Fk) = c("comcol","sameRegion","intercept")
## END Initialization ##

## BEGIN one simulation ##
source("functions/cov_TFR_fit_funcs.R")
fit_param(dimension, num_observations, matList3, sim3, adj_positions=adj_positions, Sigma=Sigma,
          filename_param_fit=paste(data_source,"sim_03_param_fit.csv",sep=""),
          filename_ests = paste(data_source,"sim_03_ests.csv",sep=""),
          filename_bic = paste(data_source,"sim_03_bic.csv",sep=""),
          filename_error_measures = paste(data_source,"sim_03_error_measures.csv",sep=""),
          compute_WSCE = TRUE)

read.csv(file=paste(data_source,"sim_03_error_measures.csv",sep=""))

# WARNING: THE MATRICES ARE NOT COMBINED THIS IS JUST A DUMMY
fit_param(dimension, num_observations, matList3, sim3, adj_positions=adj_positions, Sigma=Sigma,
          filename_param_fit=paste(data_source,"sim_03_param_fit_combined_effects.csv",sep=""),
          filename_ests = paste(data_source,"sim_03_ests_combined_effects.csv",sep=""), 
          filename_bic = paste(data_source,"sim_03_bic_combined_effects.csv",sep=""),
          filename_error_measures = paste(data_source,"sim_03_error_measures_combined_effects.csv",sep=""))
read.csv(file=paste(data_source,"sim_03_error_measures_combined_effects.csv",sep=""))

read.csv(file=paste(data_source,"sim_03_bic_combined_effects.csv",sep=""))$bic - 
  read.csv(file=paste(data_source,"sim_03_bic.csv",sep=""))$bic
## END one simulation ##

## BEGIN many simulations for mu and sigma known ##
set.seed(seed)
num_sim=40
sim_func = function(seed){
  set.seed(seed)
  sim_covs = sim_cov(num_observations, as.matrix(Sigma))
  sim_covs$dataset[is.na(dataset)] = NA
  correlation_matrix = matrix(ncol=dimension,nrow=dimension)
  
  #use pairwise correlation estimates
  for(i in (1:dimension)){
    for(j in (i:dimension)){
      dataseti_notmissing = !is.na(sim_covs$dataset[,i])
      datasetj_notmissing = !is.na(sim_covs$dataset[,j])
      correlation_matrix_ij = cor(sim_covs$dataset[dataseti_notmissing&datasetj_notmissing,i], 
                    sim_covs$dataset[dataseti_notmissing&datasetj_notmissing,j])
      correlation_matrix[i,j] = correlation_matrix_ij
      correlation_matrix[j,i] = correlation_matrix_ij
    }
  }
  sim_covs$correlation_matrix = correlation_matrix
  return(sim_covs)
} 
fit1_func = function(sim) fit_param(dimension, num_observations, matList2, sim, adj_positions=adj_positions,
                                    filename_error_measures=TRUE,save=FALSE,
                                    Sigma=Sigma, compute_WSCE = TRUE)
fit2_func = function(sim) fit_param(dimension, num_observations, matList2, sim, adj_positions=adj_positions,
                                    filename_error_measures = TRUE, save=FALSE,
                                    Sigma=Sigma, compute_WSCE = TRUE)
source("functions/cov_TFR_fit_funcs.R")
sims_errors_and_bic = sim_errors_and_bic(sim_func,fit1_func,fit2_func,num_sim)
write.csv(sims_errors_and_bic, file=paste(data_source,"sim_03_sims_errors_and_bic.csv",sep=""))
## END many simulations for mu and sigma known ##