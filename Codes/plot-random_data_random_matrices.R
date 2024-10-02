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
library(ggplot2)
library(reshape2)

data_source = "data/"
ESTS_NAMES = c("Pearson", "FM", "Glasso", "LW","IVE", "SCE", "WSCE")#c("Pearson","LW","Sparse","FM","hatSigma0","hatSigma","WSCE")
COLVEC = c("brown","grey","pink","pink3","beige","orange2","darkorange2")
PARAM1_NAMES = c("comcol", "reg", "global", "contig.beta", "contig.rho")
PARAM2_NAMES = c(PARAM1_NAMES, "comcol.and.reg", "comcol.and.config", 
                 "reg.and.config", "contig.rho")

my_theme <- theme_bw() +
  theme(strip.background = element_rect(fill = "white"), text = element_text(face="bold", size=12),
  )
theme_set(my_theme)
### Simulation 1 ###

matList2 = read_matList(filename = paste(data_source,"sim_01_matList.csv",sep=""))
(sim_01_true_param = read_param(filename=paste(data_source,"sim_01_true_param.csv",sep="")))

Sigma = CovMat_03(as.matrix(sim_01_true_param), matList2)$Sigma

# known means and variances
plot_cov(matList2,Sigma,
         SigmaHat_list=read_ests(filename=paste(data_source,"sim_01_ests.csv",sep="")),
         colvec=c("brown","grey","pink","pink3","beige","orange2","darkorange2"),
         model="corY", ests_names=ESTS_NAMES,order=c(1, 2, 3, 4, 7, 5, 6))
ggsave("atelier/sim_01_ests.jpeg", width=5.3,height=4.07,device="jpeg",dpi=700)

# unknown means and variances
plot_cov(matList2,Sigma,
         SigmaHat_list=read_ests(filename=paste(data_source,"sim_01_ests_musigma_unknown.csv",sep="")),
         colvec=c("brown","grey","pink","pink3","beige","orange2","darkorange2"),
         model="corY", ests_names=ESTS_NAMES,order=c(1, 2, 3, 4, 7, 5, 6))
ggsave("atelier/sim_01_ests_musigmaunknown.jpeg", width=5.3,height=4.07,device="jpeg",dpi=700)
# heatmap of the correlation matrix
plot_heatmaps(matList2, Sigma, 
              filename="atelier/sim_01_all_matrices_plot.jpeg")
ggsave("atelier/sim_01_all_matrices_plot.pdf", width=10.6,height=8.14)

# boxplot of error measures 
sims_errors_and_bic = read.csv(file=paste(data_source,"sim_01_sims_errors_and_bic.csv",sep=""))

param_pos = sapply(names(sims_errors_and_bic), function(s) grepl("param",s))
sims_params = sims_errors_and_bic[,param_pos]

params_1_pos = sapply(names(sims_params), function(s) grepl("param1",s))
sims_params1 = sims_params[,params_1_pos]
names(sims_params1) = PARAM1_NAMES

# errors measures of the parameters
p = 10
id_min = 1:nrow(Sigma)
plot_param_sims("atelier/sim_01_param_sims.pdf",
                sims_params1,p,Sigma,matList2,sim_01_true_param,id_min,type="Chebyshef")
# [1] "normal confidence intervals"
# comcol         reg      global contig.beta 
# 0.875       0.925       0.950       0.400 
# [1] 0.7875
# [1] "Chebyshef confidence intervals"
# comcol         reg      global contig.beta 
# 1.00        1.00        1.00        0.85 
# [1] 0.9625

sims_errors_and_bic = sims_errors_and_bic[,!param_pos]
plot_sims(sims_errors_and_bic=sims_errors_and_bic, filename="atelier/sim_01_error_measures.pdf")

# boxplot of error measures (mu and sigma unknown)
sims_errors_and_bic = read.csv(file=paste(data_source,"sim_01_sims_errors_and_bic_musigma_unknown.csv",sep=""))

param_pos = sapply(names(sims_errors_and_bic), function(s) grepl("param",s))
sims_params = sims_errors_and_bic[,param_pos]

params_1_pos = sapply(names(sims_params), function(s) grepl("param1",s))
sims_params1 = sims_params[,params_1_pos]
names(sims_params1) = PARAM1_NAMES

# errors measures of the parameters
p = 10
id_min = 1:nrow(Sigma)
plot_param_sims("atelier/sim_01_param_sims_mu_sigma_unknown.pdf",
                sims_params1,p,Sigma,matList2,sim_01_true_param,id_min,type="Chebyshef")
# [1] "normal confidence intervals"
# comcol         reg      global contig.beta 
# 0.850       0.850       0.850       0.175 
# [1] 0.68125
# [1] "Chebyshef confidence intervals"
# comcol         reg      global contig.beta 
# 0.975       0.975       0.975       0.650 
# [1] 0.89375

sims_errors_and_bic = sims_errors_and_bic[,!param_pos]
plot_sims(sims_errors_and_bic=sims_errors_and_bic, filename="atelier/sim_01_error_measures_musigma_unknown.pdf")
