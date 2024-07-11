

rm(list=ls())

#source('FunctionsCovTFR.R')
source("functions/FunctionsCovTFR_02.R")
source("functions/funcs_frob_norm_opt.R")
source("functions/cov_TFR_sim_funcs.R")
source("functions/cov_TFR_fit_funcs.R")
source("functions/cov_TFR_plot_funcs.R")
source("functions/cov_TFR_data_funcs.R")
library(corrplot)
library(viridis)

data_source = "data/"
ESTS_NAMES = c("Pearson", "FM", "Glasso", "LW","IVE", "SCE", "WSCE")#c("Pearson","LW","Sparse","FM","hatSigma0","hatSigma","WSCE")
COLVEC = c("brown","grey","pink","pink3","beige","orange2","darkorange2")
PARAM1_NAMES = c("comcol", "reg", "global", "contig.beta", "contig.rho")
PARAM2_NAMES = c(PARAM1_NAMES, "comcol.and.reg", "comcol.and.config", 
                 "reg.and.config", "contig.rho")

source("functions/get_data_covar.R")

### Read in data ###
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

parm = c(.05,0.09,.11,.74,rho)

test = calc_tilde_G_inv(matList2$Ml[[1]],matList2$Al[[1]],rho)[id_min,id_min]
diag(test)=0
max(test)*.74#maximum neighbor effect around .26
G_inv = calc_tilde_G_inv(matList2$Ml[[1]],matList2$Al[[1]],rho)[id_min,id_min]
A = matList2$Al[[1]][id_min,id_min]
A[is.na(A)] = 0 # since islands can return NA-values

covMat <- CovMat_03(parm, matList2,id_min=id_min)
Sigma <- covMat$Sigma#cov2cor(covMat$Sigma)
sim2 = sim_cov(p, as.matrix(Sigma))

### End read data ###

my_theme <- theme_bw() +
  theme(strip.background = element_rect(fill = "white"), text = element_text(face="bold", size=12),
  )
theme_set(my_theme)

# known means and variances
plot_cov(matList2,Sigma,
         SigmaHat_list=read_ests(filename=paste(data_source,"sim_03_ests.csv",sep="")),
         colvec=c("brown","grey","pink","pink3","beige","orange2","darkorange2"),
         model="corY", ests_names=ESTS_NAMES,order=c(1, 2, 3, 4, 7, 5, 6))
ggsave("atelier/sim_03_ests.jpeg", width=5.3,height=4.07,device="jpeg")


### Simulation 3 ###

matList3 = read_matList(filename = paste(data_source,"sim_02_matList.csv",sep=""))
(sim_03_true_param = read_param(filename=paste(data_source,"sim_02_true_param.csv",sep="")))
# plot_param(matList3, filename_true_param=paste(data_source,"sim_02_true_param.csv",sep=""),
#            filename_param_fit=paste(data_source,"sim_03_param_fit.csv",sep=""), id_min=id_min)

Sigma = CovMat_03(as.matrix(sim_03_true_param), matList2, id_min=id_min)$Sigma
ests = read_ests(filename=paste(data_source,"sim_03_ests.csv",sep=""))

sims_errors_and_bic = read.csv(file=paste(data_source,"sim_03_sims_errors_and_bic.csv",sep=""))

param_pos = sapply(names(sims_errors_and_bic), function(s) grepl("param",s))
sims_params = sims_errors_and_bic[,param_pos]

params_1_pos = sapply(names(sims_params), function(s) grepl("param1",s))
sims_params1 = sims_params[,params_1_pos]
names(sims_params1) = PARAM1_NAMES
library(reshape2)
ggplot(melt(sims_params1),aes(x=variable,y=value)) + geom_boxplot() + 
  geom_abline(slope=0, intercept=as.numeric(sim_03_true_param[1])) +
  geom_abline(slope=0, intercept=as.numeric(sim_03_true_param[2])) +
  geom_abline(slope=0, intercept=as.numeric(sim_03_true_param[3])) +
  geom_abline(slope=0, intercept=as.numeric(sim_03_true_param[4])) +
  geom_abline(slope=0, intercept=as.numeric(sim_03_true_param[5])) 
ggsave("atelier/sim_03_param_error_measures.pdf", width=5.3,height=4.07)

sims_errors_and_bic = sims_errors_and_bic[,!param_pos]
plot_sims(sims_errors_and_bic=sims_errors_and_bic, filename="atelier/sim_03_error_measures.pdf")

### TODO: REMOVE THE PART BELOW ###

df=read.csv(file=paste(data_source,"sim_02_modelmisspec06.csv",sep=""),row.names=1)

#df1 = df[(df$est=="hatSigma")&(df$type =="MAE"),]
df$lambda=c(" ",(df$lambda))[-1]
#df[(df$est=="Sparse"),]$value = savepoint

# renaming and changing the order
df[(df$est=="hatSigma"),]$est = "ZSCE"
df[(df$est=="hatSigma0"),]$est = "IVE"
df[(df$est=="new"),]$est = "ZZWSCE"
df[(df$est=="Sparse"),]$est = "Glasso"


df$MAE=df$value
df$estimator=df$est
plot_without_lines = ggplot(df, aes(x=lambda, y=MAE, fill=estimator)) + geom_boxplot()
for(i in (1:10)){
  plot_without_lines = plot_without_lines + geom_vline(xintercept = i+.5)
}
plot_without_lines+scale_fill_hue( labels =  expression("Glasso","IVE","LW","Pearson","SCE","WSCE"))
ggsave("atelier/sim_02_model_misspec_full.pdf", width=5.3,height=4.07,device="pdf")

# Version 2: Missing information
my_theme <- theme_bw() +
  theme(strip.background = element_rect(fill = "white"), text = element_text(face="bold", size=12),
        axis.ticks.x=element_blank()
  )
theme_set(my_theme)

df=read.csv(file=paste(data_source,"sim_02_modelmisspec_withWSCE_withmissingvalue.csv",sep=""),row.names=1)

#df1 = df[(df$est=="hatSigma")&(df$type =="MAE"),]
df$lambda=c(" ",(df$lambda))[-1]
df[(df$est=="hatSigma"),]$est = "ZSCE"
df[(df$est=="hatSigma0"),]$est = "IVE"
df[(df$est=="new"),]$est = "ZZWSCE"
#df[(df$est=="Sparse"),]$est = "Glasso"
df$MAE=df$value
df$estimator=df$est


plot_without_lines = ggplot(df, aes(x=lambda, y=MAE, fill=estimator)) + geom_boxplot()
for(i in (1:10)){
  plot_without_lines = plot_without_lines + geom_vline(xintercept = i+.5)
}
plot_without_lines + scale_fill_hue( labels = expression("Imputed","IVE","SCE","WSCE"))
ggsave("atelier/sim_02_model_misspec_partial_missing_values.pdf", width=5.3,height=4.07,device="pdf")

lambda_df=read.csv(file=paste(data_source,"sim_02_modelmisspec_withWSCE_withmissingvalue_lambdas.csv",sep=""),row.names=1)
names(lambda_df) = c("lambda","value")
lambda_df$lambda = c("",lambda_df$lambda)[-1]
ggplot(lambda_df, aes(x=lambda, y=value)) + geom_boxplot()
ggsave("atelier/sim_02_model_misspec_partial_missing_values_lambdas.pdf", width=5.3,height=4.07,device="pdf")
