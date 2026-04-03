rm(list=ls())
source("functions/FunctionsCovTFR_02.R")
source("functions/funcs_frob_norm_opt.R")
source("functions/cov_TFR_sim_funcs.R")
source("functions/cov_TFR_fit_funcs.R")
source("functions/cov_TFR_plot_funcs.R")
source("functions/cov_TFR_data_funcs.R")
data_source = "data/"

## Path 2: Use ALL countries

## BEGIN Initialization ##
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

seed=3
set.seed(seed)

dataset = FITcomps_std_total[2:12,which(all_min==1)]
dataset = dataset[,sapply(adj_positions,function(id) which(preproc_res$FITcomps_std_iso==preproc_res$iso_id_key[id]))]
dimension=dim(dataset)[2]
correlation_matrix = matrix(ncol=dimension,nrow=dimension)

#use pairwise correlation estimates
for(i in (1:dimension)){
  for(j in (i:dimension)){
    #browser()
    dataseti_notmissing = !is.na(dataset[,i])
    datasetj_notmissing = !is.na(dataset[,j])
    correlation_matrix_ij = cor(dataset[dataseti_notmissing&datasetj_notmissing,i], 
                  dataset[dataseti_notmissing&datasetj_notmissing,j])
    correlation_matrix[i,j] = correlation_matrix_ij
    correlation_matrix[j,i] = correlation_matrix_ij
  }
}
mean(is.na(correlation_matrix))
dimension = dim(matList_final$Fk[[1]])[1]
num_observations = 11
sim_final = list(correlation_matrix=correlation_matrix,dataset=dataset)

source("functions/cov_TFR_fit_funcs.R")
names(matList_final$Fk) = c("comcol","reg","intercept")
## END Initialization ##

## BEGIN check identifiability (this may take a long time) ##

# grid-search for beta/=beta' 
# (independence of the matrices is checked automatically)
beta=1.5
num_grid_sim=100
xi = (1:(num_grid_sim+1))/(num_grid_sim+1)
#tan-hyperbolic-spaced grid because beta approaches 1
beta_vec = (1-tanh(beta*(1+xi))/tanh(beta))/(min((1-tanh(beta*(1+xi))/tanh(beta))))
beta_vec = beta_vec[-length(beta_vec)]

library(lpSolve)

for(beta_1 in beta_vec){
  for(beta_2 in beta_vec){
    if(beta_1!=beta_2){
      K = length(matList_final$Fk)
      f.obj = c(rep(0,K),1,rep(0,K),1)
      
      
      
      G_matrix_1 = tilde_G_inv(M=matList_final$Ml[[1]],A=matList_final$Al[[1]],
                                    beta=beta_1)[adj_positions,adj_positions]
      G_matrix_2 = tilde_G_inv(M=matList_final$Ml[[1]],A=matList_final$Al[[1]],
                                    beta=beta_2)[adj_positions,adj_positions]
      
      for(k in 1:K){
        diag(matList_final$Fk[[k]])=1
      }
      
      # matrix combination constraints
      f.con = cbind(sapply(1:K,function(k) c(as.matrix(matList_final$Fk[[k]]))),
                    c(as.matrix(G_matrix_1)),
                    -sapply(1:K,function(k) c(as.matrix(matList_final$Fk[[k]]))),
                    -c(as.matrix(G_matrix_2)))
      f.rhs = rep(0,dim(f.con)[1])
      
      # simplex constraints
      f.con = rbind(f.con, c(rep(1,K+1),rep(0,K+1)))
      f.con = rbind(f.con, c(rep(0,K+1),rep(1,K+1)))
      f.dir = rep("=", dim(f.con)[1])
      f.rhs = c(f.rhs,1,1)
      
      # positivity constraints
      f.con = rbind(f.con, diag(2*K+2))
      f.dir = c(f.dir,rep(">=", 2*K+2))
      f.rhs = c(f.rhs, rep(0,2*K+2))
      
      solution <- lp("max",f.obj,f.con,f.dir,f.rhs)
      if(solution$objval!=0){
        print("WARNING: NOT IDENTIFIABLE")
        print(solution$solution)
      }
    }
  }
}

### END: check identifiability ###

## BEGIN fit data without interaction effects ##
fit_param(dimension, num_observations, matList_final, sim_final, adj_positions=adj_positions,
          filename_param_fit=paste(data_source,"sim_final_n195_param_fit.csv",sep=""),
          filename_ests = paste(data_source,"sim_final_n195_ests.csv",sep=""),
          filename_bic = paste(data_source,"sim_final_n195_bic.csv",sep=""),
          already_demeaned = TRUE)
read.csv(paste(data_source,"sim_final_n195_param_fit.csv",sep=""))
res_basemodel = fit_param(dimension, num_observations, matList_final, sim_final, 
                          adj_positions=adj_positions,
                          save=FALSE,
                          already_demeaned = TRUE)
write.csv(as.data.frame(res_basemodel[[2]][[6]]),
          file = "data/sim_final_n195_sigma.csv")
read.csv("data/sim_final_n195_sigma.csv")
res_basemodel[[4]]
## END fit data without interaction effects ##

## BEGIN fit data for all possible models and compute BIC ##
this_list = list()
this = c(0,0,0,0,0,0,0)

for(k in (1:(2^length(this)-1))){
  # order: comcol, reg, intercept
  pairwise_matrices_selected = (as.numeric(intToBits(k))[1:3] == 1)
  adjacency_matrix_selected = (as.numeric(intToBits(k))[4] == 1)

  # order: comcol+reg, comcol+contig, reg+contig
  combinations_selected = (as.numeric(intToBits(k))[5:7] == 1)
  
  if(combinations_selected[1]){
    if((!pairwise_matrices_selected[1])|(!pairwise_matrices_selected[2])){
      next # interaction effects can only exist if individual effects do
    }
  }
  
  if(combinations_selected[2]){
    if((!pairwise_matrices_selected[1])|(!adjacency_matrix_selected)){
      next # interaction effects can only exist if individual effects do
    }
  }
  
  if(combinations_selected[3]){
    if((!pairwise_matrices_selected[2])|(!adjacency_matrix_selected)){
      next # interaction effects can only exist if individual effects do
    }
  }
  interaction_effects_full = list(c("comcol","reg"), 
                               c("comcol","spatial"), 
                               c("reg","spatial"))
  if(sum(pairwise_matrices_selected)==0){
    matList_temp = matList_final
    matList_temp$Fk = NULL
    this0 = fit_param(dimension, num_observations, matList_temp, sim_final, 
                      adj_positions=adj_positions, save=FALSE,
                      interaction_effects = interaction_effects_full[combinations_selected],
                      already_demeaned = TRUE)
  } else if(!adjacency_matrix_selected){
    matList_temp = matList_final
    matList_temp$Fk = matList_final$Fk[pairwise_matrices_selected]
    matList_temp$vois = NULL
    this0 = fit_param(dimension, num_observations, matList_temp, sim_final, 
                      adj_positions=adj_positions, save=FALSE,
                      interaction_effects = interaction_effects_full[combinations_selected],
                      already_demeaned = TRUE)
  } else{
    matList_temp = matList_final
    matList_temp$Fk = matList_final$Fk[pairwise_matrices_selected]
    this0 = fit_param(dimension, num_observations, matList_temp, sim_final, 
                      adj_positions=adj_positions, save=FALSE,
                      interaction_effects = interaction_effects_full[combinations_selected],
                      already_demeaned = TRUE)
  }
  this_list[[k]] = this0
  print(k)
  print(as.numeric(intToBits(k))[1:7])
}
## END fit data for all possible models and compute BIC ##

## BEGIN find model with best BIC ##

# This is the bic for the model with all but the interaction effects,
# since 15 in bits is equal to 1111000
simple_model_bic = this_list[[15]][[3]]

# center the bics around the bic of the simple model
bics = sapply(this_list,function(t) t[[3]])
bics[sapply(bics,function(x) is.null(x))] = NA
bics = simplify2array(bics)
bics = bics - simple_model_bic

mat = matrix(nrow=length(bics),ncol=7)
names_of_entries = names(this_list[[length(this_list)]][[4]])
plot(sort(bics),xlab="Index of the model",ylab="BIC",main="sorted BIC values")

for(it in 1:ncol(mat)){
  vec = sapply(this_list,function(t) t[[4]][names(t[[4]])==names_of_entries[it]])
  vec[sapply(vec, function(x) is.null(x) | (length(x)==0))] = NA
  mat[,it] = simplify2array(vec)
}
colnames(mat) = c(paste0(1:7,names_of_entries))
rownames(mat) = round(bics,1)
mat = mat[rowSums(is.na(mat))<ncol(mat),]
df = as.data.frame(mat)
write.csv(df,file=paste(data_source,"sim_final_n195_model_choice.csv",sep=""))

ix_optim_model = which.min(bics)
model_params=as.numeric(intToBits(ix_optim_model))[1:7]
model_params # everything but "comcol and reg"

fit_param(dimension, num_observations, matList_final, sim_final, adj_positions=adj_positions,
          filename_param_fit=paste(data_source,"sim_final_n195_combined_param_fit.csv",sep=""),
          filename_ests = paste(data_source,"sim_final_n195_combined_ests.csv",sep=""),
          filename_bic = paste(data_source,"sim_final_n195_combined_bic.csv",sep=""),
          interaction_effects = list(c("comcol","spatial"), 
                                     c("reg","spatial")),
          compute_WSCE=TRUE,
          already_demeaned = TRUE)
read.csv(paste(data_source,"sim_final_n195_combined_param_fit.csv",sep=""))
# the WSCE was equal to the SCE
fit_param(dimension, num_observations, matList_final, sim_final, adj_positions=adj_positions,
          filename_param_fit=paste(data_source,"sim_final_n195_combined_param_fit.csv",sep=""),
          filename_ests = paste(data_source,"sim_final_n195_combined_ests.csv",sep=""),
          filename_bic = paste(data_source,"sim_final_n195_combined_bic.csv",sep=""),
          interaction_effects = list(c("comcol","spatial"), 
                                     c("reg","spatial")),
          already_demeaned = TRUE)
parm=read.csv(paste(data_source,"sim_final_n195_combined_param_fit.csv",sep=""))

res = fit_param(dimension, num_observations, matList_final, sim_final, adj_positions=adj_positions,
                interaction_effects = list(c("comcol","spatial"), 
                                           c("reg","spatial")),
                save = FALSE,
                already_demeaned = TRUE)
write.csv(as.data.frame(res[[2]][[6]]),
          file = "data/sim_final_n195_sigma_with_interactions.csv")
res[[4]]
## END find model with best BIC ##