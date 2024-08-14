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
id_min = preproc_res$id_min

# Parms
n <- 195; p <- 11; rho = .35

matList2 = matList_final
write_matList(matList=matList2, filename=paste(data_source,"sim_02_matList.csv",sep=""))

parm = c(.05,0.09,.11,.74,rho)
write_param(t(parm),matList=matList2,
            filename=paste(data_source,"sim_02_true_param.csv",sep=""))

test = calc_tilde_G_inv(matList2$Ml[[1]],matList2$Al[[1]],rho)[id_min,id_min]
diag(test)=0
max(test)*.74#maximum neighbor effect around .26
G_inv = calc_tilde_G_inv(matList2$Ml[[1]],matList2$Al[[1]],rho)[id_min,id_min]
A = matList2$Al[[1]][id_min,id_min]
A[is.na(A)] = 0 # since islands can return NA-values

covMat <- CovMat_03(parm, matList2,id_min=id_min)
Sigma <- covMat$Sigma#cov2cor(covMat$Sigma)
sim2 = sim_cov(p, as.matrix(Sigma))

source("functions/cov_TFR_fit_funcs.R")
fit_param(n, p, matList2, sim2, Sigma=Sigma, id_min=id_min,
          filename_param_fit=paste(data_source,"sim_02_param_fit.csv",sep=""),
          filename_ests = paste(data_source,"sim_02_ests.csv",sep=""),
          filename_bic = paste(data_source,"sim_02_bic.csv",sep=""),
          filename_error_measures = paste(data_source,"sim_02_error_measures.csv",sep=""),
          compute_WSCE = TRUE)
read.csv(file=paste(data_source,"sim_02_error_measures.csv",sep=""))

fit_param(n, p, matList2, sim2, Sigma=Sigma, id_min=id_min,
          filename_param_fit=paste(data_source,"sim_02_param_fit_musigma_unknown.csv",sep=""),
          filename_ests = paste(data_source,"sim_02_ests_musigma_unknown.csv",sep=""),
          filename_bic = paste(data_source,"sim_02_bic_musigma_unknown.csv",sep=""),
          filename_error_measures = paste(data_source,"sim_02_error_measures_musigma_unknown.csv",sep=""),
          compute_WSCE = TRUE, normalize_data = TRUE)
read.csv(file=paste(data_source,"sim_02_error_measures_musigma_unknown.csv",sep=""))

set.seed(seed)
num_sim=10
sim_func = function(seed) sim_cov_with_seed(p, as.matrix(Sigma), seed)
fit1_func = function(sim) fit_param(n, p, matList2, sim, id_min=id_min,
                                    filename_error_measures=TRUE,save=FALSE,Sigma=as.matrix(Sigma),
                                    compute_WSCE = TRUE)
fit2_func = function(sim) fit_param(n, p, matList2, sim, id_min=id_min,
                                    filename_error_measures = TRUE,
                                    link=combined_matList, link_der_rho = link_der_combined,
                                    save=FALSE,Sigma=as.matrix(Sigma), compute_WSCE = TRUE)
sims_errors_and_bic = sim_errors_and_bic(sim_func,fit1_func,fit2_func,num_sim)
write.csv(sims_errors_and_bic,
          file=paste(data_source,"sim_02_sims_errors_and_bic.csv",sep=""))

set.seed(seed)
num_sim=10
sim_func = function(seed) sim_cov_with_seed(p, as.matrix(Sigma), seed)
fit1_func = function(sim) fit_param(n, p, matList2, sim, id_min=id_min,
                                    filename_error_measures=TRUE,save=FALSE,Sigma=as.matrix(Sigma),
                                    compute_WSCE = TRUE,normalize_data = TRUE)
fit2_func = function(sim) fit_param(n, p, matList2, sim, id_min=id_min,
                                    filename_error_measures = TRUE,
                                    link=combined_matList, link_der_rho = link_der_combined,
                                    save=FALSE,Sigma=as.matrix(Sigma), compute_WSCE = TRUE, normalize_data = TRUE)
sims_errors_and_bic = sim_errors_and_bic(sim_func,fit1_func,fit2_func,num_sim)
write.csv(sims_errors_and_bic,
          file=paste(data_source,"sim_02_sims_errors_and_bic_musigma_unknown.csv",sep=""))

## Compute performance with varying n

# read in data about regions
library(cepiigeodist)

covar_ord_by_id <- covar %>% arrange(id_col) %>% group_by(id_col) %>% 
  slice(1) %>%ungroup()
all(covar_ord_by_id%>%  pull(name_o) == names_by_id) # on retrouve bien comme Martin:

iso3 <- covar_ord_by_id %>% pull(code_o)
missing_iso <- iso3[!iso3 %in%geo_cepii$iso3 ]

## Let's get the information (of  colonizer) using covar 

### Récupération des informations dans geo_cepii
get_count_info <- function(iso_sel, missing_iso, covar, geo_cepii){
  
  same_count <-  covar  %>% filter(code_o ==iso_sel &  !code_d%in% missing_iso) %>% 
    filter(comcol==1)  %>% slice(1) %>%  arrange(id_col) %>% pull(code_d)
  
  if(length(same_count)!=0){
    tibble(iso = iso_sel, iso3 = same_count) %>%
      left_join(geo_cepii %>% select(iso3, colonizer1)) %>% 
      unique() %>% pull(colonizer1)
    
  }else{
    "."
  }
  
}

colon_full<- tibble(iso3 = missing_iso) %>% 
  mutate(colonizer1 = map_chr(iso3, ~get_count_info(., missing_iso, covar, geo_cepii))) %>% 
  bind_rows(geo_cepii %>% select(iso3,colonizer1))


## Now let's merge everything we need 

info_coutries <- countries_num %>% 
  select(country_name,`ISO3 Alpha-code`, reg_name) %>%
  left_join(colon_full, by =c("ISO3 Alpha-code" = "iso3"))  %>% 
  unique()

Sigma = as.matrix(Sigma)
Sigma_countries <- Sigma %>% 
  as.data.frame() %>% 
  mutate(iso3 =  iso3[id_min] ) %>% 
  left_join(info_coutries, by =c( "iso3" = "ISO3 Alpha-code"))
Sigma_countries$reg_name

# end read in data about regions

#Check number of countries per index set
index_14 = which((Sigma_countries$reg_name=="Middle Africa")|(Sigma_countries$reg_name=="Southern Africa"))
index_32 = which((Sigma_countries$reg_name=="Eastern Africa")|(Sigma_countries$reg_name=="Southern Africa")|(Sigma_countries$reg_name=="Middle Africa"))
index_55 = which((Sigma_countries$reg_name=="Eastern Africa")|(Sigma_countries$reg_name=="Middle Africa")|
      (Sigma_countries$reg_name=="Western Africa")|(Sigma_countries$reg_name=="Southern Africa")|
      (Sigma_countries$reg_name=="Northern Africa")|(Sigma_countries$reg_name=="Northern America")|
        (Sigma_countries$reg_name=="Central America")|(Sigma_countries$reg_name=="Southern America"))
index_103 = which((Sigma_countries$reg_name=="Eastern Africa")|(Sigma_countries$reg_name=="Middle Africa")|
      (Sigma_countries$reg_name=="Western Africa")|(Sigma_countries$reg_name=="Southern Africa")|
      (Sigma_countries$reg_name=="Northern Africa")|
      (Sigma_countries$reg_name=="Eastern Asia") | (Sigma_countries$reg_name=="South-Central Asia") |
      (Sigma_countries$reg_name=="South-Eastern Asia")|(Sigma_countries$reg_name=="Western Asia")|(Sigma_countries$reg_name=="Northern America")|
        (Sigma_countries$reg_name=="Central America")|(Sigma_countries$reg_name=="Southern America"))

index_n = list(index_14,index_32,index_55,index_103,1:195)

set.seed(seed)
res_list=list()
matList1=matList2
source("functions/cov_TFR_fit_funcs.R")

res_list_list = list()
lambdas = matrix(ncol=length(index_n),nrow=40)

for(N in 1:40){
  for(i in seq_along(index_n)){
    
    Sigma_stretched = Sigma[index_n[[i]],index_n[[i]]]
    sim2 = sim_cov_with_seed(p, as.matrix(Sigma),seed=N)
    sim_stretched = sim2
    sim_stretched$Y = sim2$Y[,index_n[[i]]]
    sim_stretched$corY = cor_from_standard_errors(sim_stretched$Y)

    for(k in (1:3)){
      matList1$Fk[[k]] = matList2$Fk[[k]][index_n[[i]],index_n[[i]]]
    }
    res_list[[i]] = c(fit_param(length(index_n[[i]]), p, matList1, sim_stretched, 
                                id_min=id_min[index_n[[i]]], save=FALSE,
                                filename_error_measures=TRUE,
                                Sigma=Sigma_stretched, compute_WSCE = TRUE),list(sim_stretched))
    #browser()
    
  }
  res_list_list[[N]] = res_list
}

# Only save the error measures, lambda, and the parameters.
read_data_varying_n = function(N){
  t(sapply(seq_along(index_n), function(i) c(res_list_list[[N]][[i]][[1]][1,],
                                             res_list_list[[N]][[i]][[4]],
                                             res_list_list[[N]][[i]][[5]],
                                             paste0(length(index_n[[i]])))))
}

data=cbind(read_data_varying_n(1),1)
for(i in (2:N)){
  data = rbind(data, cbind(read_data_varying_n(i),i))
}
write.csv(as.data.frame(data),file=paste(data_source,"sim_02_different_n.csv",sep=""))

### model misspecification ###

matList3 = read_matList(filename = paste(data_source,"sim_01_matList.csv",sep=""))

set.seed(seed)

# parameter for Bar(Sigma) (rho does not matter since spatial effect is 0)
(sim_01_true_param = c(.99,0,0,0,.982))#.982
Sigma_extra = CovMat_03(as.matrix(sim_01_true_param), matList3)$Sigma

set.seed(seed)
Y = FITcomps_std_total[2:12,which(all_min==1)]
Y = Y[,sapply(id_min,function(id) which(preproc_res$FITcomps_std_iso==preproc_res$iso_id_key[id]))]

# Compute performance over model misspecification
extra_matrix = Sigma_extra[1:195,1:195]#cov2cor(rWishart(1,n,Sigma=Sigma2)[,,1])
mean(abs(extra_matrix[extra_matrix!=1]))
lambda_vec = (0:10)/10
res_list=list()
matList1=matList_final
source("functions/cov_TFR_fit_funcs.R")

res_list_list = list()

est_new = list()
est_new_res = list()
N_samples = 10
lambdas = matrix(ncol=length(lambda_vec),nrow=N_samples)

for(N in 1:N_samples){
  for(i in seq_along(lambda_vec)){
    Sigma_stretched = lambda_vec[i]*extra_matrix + (1-lambda_vec[i])*as.matrix(Sigma)
    sim_stretched = sim_cov_with_seed(p, Sigma_stretched,seed=N)
    
    res_list[[i]] = c(fit_param(n, p, matList1, sim_stretched, 
                                id_min=id_min, save=FALSE,
                                filename_error_measures=TRUE,
                                Sigma=Sigma_stretched,compute_WSCE = TRUE),list(sim_stretched))
    lambdas[N,i]=res_list[[i]][[4]]
    
    est_new[[i]] = 0
    est_new_res[[i]] = c(mean(abs(est_new[[i]]- Sigma_stretched)),sqrt(mean((est_new[[i]]- Sigma_stretched)^2)))
    res_list[[i]][[1]] = cbind(res_list[[i]][[1]],est_new_res[[i]])
    res_list[[i]][[3]][[length(res_list[[i]][[3]]) + 1]] = est_new[[i]]
  }
  res_list_list[[N]] = res_list
}
list_list = list(lambdas,res_list_list)

data = t(sapply(1:length(res_list_list),
                function(t) sapply(1:length(res_list_list[[1]]), 
                                   function(s) res_list_list[[t]][[s]][[1]][1,1])))
for(i in 1:2){
  for(j in 1:8){
    if(mean(c(i,j)==c(1,1))<1){
      data=cbind(data,t(sapply(1:length(res_list_list),
                               function(t) sapply(1:length(res_list_list[[1]]), 
                                                  function(s) res_list_list[[t]][[s]][[1]][i,j]))))
    }
  }
}

names_ests = c("Pearson","LW","Sparse","FM","BestWSCE","hatSigma0","hatSigma","new")
df=as.data.frame(data)
names(df)=c(paste("MAE",c(sapply(1:8,function(s) paste(lambda_vec,rep(names_ests[s],length(lambda_vec)))))),
            paste("RMSE",c(sapply(1:8,function(s) paste(lambda_vec,rep(names_ests[s],length(lambda_vec)))))))
library(reshape2)
df=melt(df)
df$est = apply(as.matrix(df$variable),1,function(s) strsplit(s," ")[[1]][3])
df$lambda = apply(as.matrix(df$variable),1,function(s) strsplit(s," ")[[1]][2])
df$type = apply(as.matrix(df$variable),1,function(s) strsplit(s," ")[[1]][1])
df$variable=NULL
write.csv(df,file=paste(data_source,"sim_02_modelmisspec_withWSCE.csv",sep=""))
write.csv(lambdas,file=paste(data_source,"sim_02_modelmisspec_withWSCE_lambdas.csv",sep=""))

### model misspecification (missing errors) ###
matList3 = read_matList(filename = paste(data_source,"sim_01_matList.csv",sep=""))

set.seed(seed)
(sim_01_true_param = c(.99,0,0,0,.982)) # same as above
Sigma_extra = CovMat_03(as.matrix(sim_01_true_param), matList3)$Sigma

set.seed(seed)
Y = FITcomps_std_total[2:12,which(all_min==1)]
Y = Y[,sapply(id_min,function(id) which(preproc_res$FITcomps_std_iso==preproc_res$iso_id_key[id]))]

# Compute performance over model misspecification
extra_matrix = Sigma_extra[1:195,1:195]
mean(abs(extra_matrix[extra_matrix!=1]))
lambda_vec = (0:10)/10
res_list=list()
matList1=matList_final

source("functions/cov_TFR_fit_funcs.R")

res_list_list = list()

est_new = list()
est_new_res = list()
N_samples = 10
lambdas = matrix(ncol=length(lambda_vec),nrow=N_samples)

for(N in 1:N_samples){
  for(i in seq_along(lambda_vec)){
    Sigma_stretched = lambda_vec[i]*extra_matrix + (1-lambda_vec[i])*as.matrix(Sigma)
    sim_stretched = sim_cov_with_seed(p, Sigma_stretched,seed=N)
    sim_stretched$Y[is.na(Y)] = NA

    res_list[[i]] = c(fit_param(n, p, matList1, sim_stretched, 
                                id_min=id_min, save=FALSE,
                                filename_error_measures=TRUE,
                                Sigma=Sigma_stretched,compute_WSCE = TRUE),list(sim_stretched))
    lambdas[N,i]=res_list[[i]][[4]]
    
    est_new[[i]] = 0  
    est_new_res[[i]] = c(mean(abs(est_new[[i]]- Sigma_stretched)),sqrt(mean((est_new[[i]]- Sigma_stretched)^2)))
    res_list[[i]][[1]] = cbind(res_list[[i]][[1]],est_new_res[[i]])
    res_list[[i]][[3]][[length(res_list[[i]][[3]]) + 1]] = est_new[[i]]
  }
  res_list_list[[N]] = res_list
}
list_list = list(lambdas,res_list_list)

data = t(sapply(1:length(res_list_list),
                function(t) sapply(1:length(res_list_list[[1]]), 
                                   function(s) res_list_list[[t]][[s]][[1]][1,1])))
for(i in 1:2){
  for(j in 1:8){
    if(mean(c(i,j)==c(1,1))<1){
      data=cbind(data,t(sapply(1:length(res_list_list),
                               function(t) sapply(1:length(res_list_list[[1]]), 
                                                  function(s) res_list_list[[t]][[s]][[1]][i,j]))))
    }
  }
}

names_ests = c("Pearson","LW","Sparse","FM","BestWSCE","hatSigma0","hatSigma","new")
df=as.data.frame(data)
names(df)=c(paste("MAE",c(sapply(1:8,function(s) paste(lambda_vec,rep(names_ests[s],length(lambda_vec)))))),
            paste("RMSE",c(sapply(1:8,function(s) paste(lambda_vec,rep(names_ests[s],length(lambda_vec)))))))
library(reshape2)
df=melt(df)
df$est = apply(as.matrix(df$variable),1,function(s) strsplit(s," ")[[1]][3])
df$lambda = apply(as.matrix(df$variable),1,function(s) strsplit(s," ")[[1]][2])
df$type = apply(as.matrix(df$variable),1,function(s) strsplit(s," ")[[1]][1])
df$variable=NULL
write.csv(df,file=paste(data_source,"sim_02_modelmisspec_withWSCE_missing_values.csv",sep=""))
write.csv(lambdas,file=paste(data_source,"sim_02_modelmisspec_withWSCE_missing_values_lambdas.csv",sep=""))
