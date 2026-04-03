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

## BEGIN Initialization ##

seed=3
set.seed(seed)

source("functions/get_data_covar.R")

read_plot = read_plot_FITcomps_std(filename="../Data/TFR_pieces_202311/standardized_residuals_202311/FITcomps_std_residuals_sample%i_202311.txt") # Initializing and plotting the real data

#results of print() output:
# [1] "1958 121"
# [1] "1963 98"
# [1] "1968 74"
# [1] "1973 56"
# [1] "1978 35"
# [1] "1983 26"
# [1] "1988 7"
# [1] "1993 4"
# [1] "1998 2"
# [1] "2003 0"
# [1] "2008 0"
# [1] "percentage missing values:"
# [1] 0.1961967

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
parm = c(.05,0.09,.11,.74,beta)
# write_param(t(parm),matList=matList2,
#             filename=paste(data_source,"sim_02_true_param.csv",sep=""))

test = tilde_G_inv(matList2$Ml[[1]],matList2$Al[[1]],beta)[adj_positions,adj_positions]
diag(test)=0
max(test)*.74#maximum neighbor effect around .26
G_inv = tilde_G_inv(matList2$Ml[[1]],matList2$Al[[1]],beta)[adj_positions,adj_positions]
A = matList2$Al[[1]][adj_positions,adj_positions]
A[is.na(A)] = 0 # since islands can return NA-values
matList2$Gl[[1]] = G_inv

covMat <- CovMat_03(parm, matList2,adj_positions=adj_positions)
Sigma <- covMat$Sigma#cov2cor(covMat$Sigma)

## END Initialization ##

## BEGIN one simulation ##

sim2 = sim_cov(num_observations, as.matrix(Sigma),degrees_of_freedom = 10)
names(matList2$Fk) = c("comcol","sameRegion","intercept")
source("functions/cov_TFR_fit_funcs.R")
fit_param(dimension, num_observations, matList2, sim2, Sigma=Sigma, adj_positions=adj_positions,
          filename_param_fit=paste(data_source,"sim_04_param_fit.csv",sep=""),
          filename_ests = paste(data_source,"sim_04_ests.csv",sep=""),
          filename_bic = paste(data_source,"sim_04_bic.csv",sep=""),
          filename_error_measures = paste(data_source,"sim_04_error_measures.csv",sep=""),
          compute_WSCE = TRUE)
read.csv(file=paste(data_source,"sim_04_error_measures.csv",sep=""))

fit_param(dimension, num_observations, matList2, sim2, Sigma=Sigma, adj_positions=adj_positions,
          filename_param_fit=paste(data_source,"sim_04_param_fit_musigma_unknown.csv",sep=""),
          filename_ests = paste(data_source,"sim_04_ests_musigma_unknown.csv",sep=""),
          filename_bic = paste(data_source,"sim_04_bic_musigma_unknown.csv",sep=""),
          filename_error_measures = paste(data_source,"sim_04_error_measures_musigma_unknown.csv",sep=""),
          compute_WSCE = TRUE, normalize_data = TRUE)
read.csv(file=paste(data_source,"sim_04_error_measures_musigma_unknown.csv",sep=""))

## END one simulation ##

## BEGIN many simulations for mu and sigma known, df=200 ##
set.seed(seed)
num_sim=40
names(matList2$Fk) = c("comcol","sameRegion","intercept")
sim_func = function(seed) sim_cov_with_seed(num_observations, as.matrix(Sigma), seed, degrees_of_freedom = 200)
fit1_func = function(sim) fit_param(dimension, num_observations, matList2, sim, adj_positions=adj_positions,
                                    filename_error_measures=TRUE,save=FALSE,Sigma=as.matrix(Sigma),
                                    compute_WSCE = TRUE)
fit2_func = function(sim) fit_param(dimension, num_observations, matList2, sim, adj_positions=adj_positions,
                                    filename_error_measures = TRUE,
                                    interaction_effects=list(c("comcol","sameRegion"),
                                                             c("comcol","spatial"),
                                                             c("sameRegion","spatial")),
                                    save=FALSE,Sigma=as.matrix(Sigma), compute_WSCE = TRUE)
sims_errors_and_bic = sim_errors_and_bic(sim_func,fit1_func,fit2_func,num_sim)
write.csv(sims_errors_and_bic,
          file=paste(data_source,"sim_04_sims_errors_and_bic_df200.csv",sep=""))
## END many simulations for mu and sigma known, df=200 ##

## BEGIN many simulations for mu and sigma known, df=500 ##
set.seed(seed)
num_sim=40
names(matList2$Fk) = c("comcol","sameRegion","intercept")
sim_func = function(seed) sim_cov_with_seed(num_observations, as.matrix(Sigma), seed, degrees_of_freedom = 500)
fit1_func = function(sim) fit_param(dimension, num_observations, matList2, sim, adj_positions=adj_positions,
                                    filename_error_measures=TRUE,save=FALSE,Sigma=as.matrix(Sigma),
                                    compute_WSCE = TRUE)
fit2_func = function(sim) fit_param(dimension, num_observations, matList2, sim, adj_positions=adj_positions,
                                    filename_error_measures = TRUE,
                                    interaction_effects=list(c("comcol","sameRegion"),
                                                             c("comcol","spatial"),
                                                             c("sameRegion","spatial")),
                                    save=FALSE,Sigma=as.matrix(Sigma), compute_WSCE = TRUE)
sims_errors_and_bic = sim_errors_and_bic(sim_func,fit1_func,fit2_func,num_sim)
write.csv(sims_errors_and_bic,
          file=paste(data_source,"sim_04_sims_errors_and_bic_df500.csv",sep=""))
## END many simulations for mu and sigma known, df=500 ##

## BEGIN many simulations for mu and sigma known, df=1000 ##
set.seed(seed)
num_sim=40
names(matList2$Fk) = c("comcol","sameRegion","intercept")
sim_func = function(seed) sim_cov_with_seed(num_observations, as.matrix(Sigma), seed, degrees_of_freedom = 1000)
fit1_func = function(sim) fit_param(dimension, num_observations, matList2, sim, adj_positions=adj_positions,
                                    filename_error_measures=TRUE,save=FALSE,Sigma=as.matrix(Sigma),
                                    compute_WSCE = TRUE)
fit2_func = function(sim) fit_param(dimension, num_observations, matList2, sim, adj_positions=adj_positions,
                                    filename_error_measures = TRUE,
                                    interaction_effects=list(c("comcol","sameRegion"),
                                                             c("comcol","spatial"),
                                                             c("sameRegion","spatial")),
                                    save=FALSE,Sigma=as.matrix(Sigma), compute_WSCE = TRUE)
sims_errors_and_bic = sim_errors_and_bic(sim_func,fit1_func,fit2_func,num_sim)
write.csv(sims_errors_and_bic,
          file=paste(data_source,"sim_04_sims_errors_and_bic_df1000.csv",sep=""))
## END many simulations for mu and sigma known, df=1000 ##

## BEGIN many simulations for mu and sigma known, df=10 ##
set.seed(seed)
num_sim=40
names(matList2$Fk) = c("comcol","sameRegion","intercept")
sim_func = function(seed) sim_cov_with_seed(num_observations, as.matrix(Sigma), seed, degrees_of_freedom = 10)
fit1_func = function(sim) fit_param(dimension, num_observations, matList2, sim, adj_positions=adj_positions,
                                    filename_error_measures=TRUE,save=FALSE,Sigma=as.matrix(Sigma),
                                    compute_WSCE = TRUE)
fit2_func = function(sim) fit_param(dimension, num_observations, matList2, sim, adj_positions=adj_positions,
                                    filename_error_measures = TRUE,
                                    interaction_effects=list(c("comcol","sameRegion"),
                                                             c("comcol","spatial"),
                                                             c("sameRegion","spatial")),
                                    save=FALSE,Sigma=as.matrix(Sigma), compute_WSCE = TRUE)
sims_errors_and_bic = sim_errors_and_bic(sim_func,fit1_func,fit2_func,num_sim)
write.csv(sims_errors_and_bic,
          file=paste(data_source,"sim_04_sims_errors_and_bic_df10.csv",sep=""))
## END many simulations for mu and sigma known, df=10 ##

## BEGIN many simulations for mu and sigma known, df=20 ##
set.seed(seed)
num_sim=40
names(matList2$Fk) = c("comcol","sameRegion","intercept")
sim_func = function(seed) sim_cov_with_seed(num_observations, as.matrix(Sigma), seed, degrees_of_freedom = 20)
fit1_func = function(sim) fit_param(dimension, num_observations, matList2, sim, adj_positions=adj_positions,
                                    filename_error_measures=TRUE,save=FALSE,Sigma=as.matrix(Sigma),
                                    compute_WSCE = TRUE)
fit2_func = function(sim) fit_param(dimension, num_observations, matList2, sim, adj_positions=adj_positions,
                                    filename_error_measures = TRUE,
                                    interaction_effects=list(c("comcol","sameRegion"),
                                                             c("comcol","spatial"),
                                                             c("sameRegion","spatial")),
                                    save=FALSE,Sigma=as.matrix(Sigma), compute_WSCE = TRUE)
sims_errors_and_bic = sim_errors_and_bic(sim_func,fit1_func,fit2_func,num_sim)
write.csv(sims_errors_and_bic,
          file=paste(data_source,"sim_04_sims_errors_and_bic_df20.csv",sep=""))
## END many simulations for mu and sigma known, df=20 ##

## BEGIN many simulations for mu and sigma known, df=50 ##
set.seed(seed)
num_sim=40
names(matList2$Fk) = c("comcol","sameRegion","intercept")
sim_func = function(seed) sim_cov_with_seed(num_observations, as.matrix(Sigma), seed, degrees_of_freedom = 50)
fit1_func = function(sim) fit_param(dimension, num_observations, matList2, sim, adj_positions=adj_positions,
                                    filename_error_measures=TRUE,save=FALSE,Sigma=as.matrix(Sigma),
                                    compute_WSCE = TRUE)
fit2_func = function(sim) fit_param(dimension, num_observations, matList2, sim, adj_positions=adj_positions,
                                    filename_error_measures = TRUE,
                                    interaction_effects=list(c("comcol","sameRegion"),
                                                             c("comcol","spatial"),
                                                             c("sameRegion","spatial")),
                                    save=FALSE,Sigma=as.matrix(Sigma), compute_WSCE = TRUE)
sims_errors_and_bic = sim_errors_and_bic(sim_func,fit1_func,fit2_func,num_sim)
write.csv(sims_errors_and_bic,
          file=paste(data_source,"sim_04_sims_errors_and_bic_df50.csv",sep=""))
## END many simulations for mu and sigma known, df=50 ##

## BEGIN many simulations for mu and sigma known, df=100 ##
set.seed(seed)
num_sim=40
names(matList2$Fk) = c("comcol","sameRegion","intercept")
sim_func = function(seed) sim_cov_with_seed(num_observations, as.matrix(Sigma), seed, degrees_of_freedom = 100)
fit1_func = function(sim) fit_param(dimension, num_observations, matList2, sim, adj_positions=adj_positions,
                                    filename_error_measures=TRUE,save=FALSE,Sigma=as.matrix(Sigma),
                                    compute_WSCE = TRUE)
fit2_func = function(sim) fit_param(dimension, num_observations, matList2, sim, adj_positions=adj_positions,
                                    filename_error_measures = TRUE,
                                    interaction_effects=list(c("comcol","sameRegion"),
                                                             c("comcol","spatial"),
                                                             c("sameRegion","spatial")),
                                    save=FALSE,Sigma=as.matrix(Sigma), compute_WSCE = TRUE)
sims_errors_and_bic = sim_errors_and_bic(sim_func,fit1_func,fit2_func,num_sim)
write.csv(sims_errors_and_bic,
          file=paste(data_source,"sim_04_sims_errors_and_bic_df100.csv",sep=""))
## END many simulations for mu and sigma known, df=100 ##
