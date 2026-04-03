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

### Simulation 2 ###

## BEGIN reading in the data ##

source("functions/get_data_covar.R")

read_plot = read_plot_FITcomps_std(filename="../Data/TFR_pieces_202311/standardized_residuals_202311/FITcomps_std_residuals_sample%i_202311.txt") # Initializing and plotting the real data
FITcomps_std_total = read_plot$FITcomps_std_total
FITcomps_std = read_plot$FITcomps_std

read_names = read_names_FITcomps_std_total(FITcomps_std_total, covar, model="all_values")
names_by_id = read_names$names_by_id
all_min = read_names$all_min

preproc_res = preproc_FITcomps_std(all_min, names_by_id, FITcomps_std, covar)
matList_final = preproc_res$matList_final
adj_positions = preproc_res$adj_positions

## END reading the data ##

# ## BEGIN plot 1 simulation ##
# matList2 = read_matList(filename = paste(data_source,"sim_02_matList.csv",sep=""))
# (sim_02_true_param = read_param(filename=paste(data_source,"sim_02_true_param.csv",sep="")))
# 
# Sigma = CovMat_03(as.matrix(sim_02_true_param), matList2, adj_positions=adj_positions)$Sigma
# ests = read_ests(filename=paste(data_source,"sim_04_ests.csv",sep=""))
# 
# # known means and variances
# plot_cov(matList2,Sigma,
#          SigmaHat_list=read_ests(filename=paste(data_source,"sim_04_ests.csv",sep="")),
#          colvec=c("brown","grey","pink","pink3","beige","orange2","darkorange2"),
#          model="correlation_matrix", ests_names=ESTS_NAMES,order=c(1, 2, 3, 4, 7, 5, 6))
# ggsave("atelier/sim_04_ests.jpeg", width=5.3,height=4.07,device="jpeg",dpi=700)
# 
# # unknown means and variances
# plot_cov(matList2,Sigma,
#          SigmaHat_list=read_ests(filename=paste(data_source,"sim_04_ests_musigma_unknown.csv",sep="")),
#          colvec=c("brown","grey","pink","pink3","beige","orange2","darkorange2"),
#          model="correlation_matrix", ests_names=ESTS_NAMES,order=c(1, 2, 3, 4, 7, 5, 6))
# ggsave("atelier/sim_04_ests_musigmaunknown.jpeg", width=5.3,height=4.07,device="jpeg")
# 
# # heatmap of the correlation matrix
# plot_heatmaps(matList2, Sigma,
#               filename=paste(data_source,"sim_02_all_matrices_plot.jpeg",sep=""))
# ggsave("atelier/sim_02_all_matrices_plot.pdf", width=10.6,height=8.14)
# ## END plot 1 simulation ##
# 
# ## BEGIN boxplots and confidence intervals ##
# num_observations <- 11
# 
# error measures
sims_errors_and_bic_df1000 = read.csv(file=paste(data_source,"sim_04_sims_errors_and_bic_df1000.csv",sep=""))
param_pos = sapply(names(sims_errors_and_bic_df1000), function(s) grepl("param",s))
sims_params = sims_errors_and_bic_df1000[,param_pos]
sims_errors_and_bic_df1000 = sims_errors_and_bic_df1000[,!param_pos]
sims_errors_and_bic_df1000$X = 1000
sims_errors_and_bic_df500 = read.csv(file=paste(data_source,"sim_04_sims_errors_and_bic_df500.csv",sep=""))
param_pos = sapply(names(sims_errors_and_bic_df500), function(s) grepl("param",s))
sims_params = sims_errors_and_bic_df500[,param_pos]
sims_errors_and_bic_df500 = sims_errors_and_bic_df500[,!param_pos]
sims_errors_and_bic_df500$X = 500
sims_errors_and_bic_df200 = read.csv(file=paste(data_source,"sim_04_sims_errors_and_bic_df200.csv",sep=""))
param_pos = sapply(names(sims_errors_and_bic_df200), function(s) grepl("param",s))
sims_params = sims_errors_and_bic_df200[,param_pos]
sims_errors_and_bic_df200 = sims_errors_and_bic_df200[,!param_pos]
sims_errors_and_bic_df200$X = 200
sims_errors_and_bic_df100 = read.csv(file=paste(data_source,"sim_04_sims_errors_and_bic_df100.csv",sep=""))
param_pos = sapply(names(sims_errors_and_bic_df100), function(s) grepl("param",s))
sims_params = sims_errors_and_bic_df100[,param_pos]
sims_errors_and_bic_df100 = sims_errors_and_bic_df100[,!param_pos]
sims_errors_and_bic_df100$X = 100
sims_errors_and_bic_df50 = read.csv(file=paste(data_source,"sim_04_sims_errors_and_bic_df50.csv",sep=""))
param_pos = sapply(names(sims_errors_and_bic_df50), function(s) grepl("param",s))
sims_params = sims_errors_and_bic_df50[,param_pos]
sims_errors_and_bic_df50 = sims_errors_and_bic_df50[,!param_pos]
sims_errors_and_bic_df50$X = 50
sims_errors_and_bic_df20 = read.csv(file=paste(data_source,"sim_04_sims_errors_and_bic_df20.csv",sep=""))
param_pos = sapply(names(sims_errors_and_bic_df20), function(s) grepl("param",s))
sims_params = sims_errors_and_bic_df20[,param_pos]
sims_errors_and_bic_df20 = sims_errors_and_bic_df20[,!param_pos]
sims_errors_and_bic_df20$X = 20
sims_errors_and_bic_df10 = read.csv(file=paste(data_source,"sim_04_sims_errors_and_bic_df10.csv",sep=""))
param_pos = sapply(names(sims_errors_and_bic_df10), function(s) grepl("param",s))
sims_params = sims_errors_and_bic_df10[,param_pos]
sims_errors_and_bic_df10 = sims_errors_and_bic_df10[,!param_pos]
sims_errors_and_bic_df10$X = 10
sims_errors_and_bic_dffull = rbind(sims_errors_and_bic_df10,
                                   sims_errors_and_bic_df20,
                                   sims_errors_and_bic_df50,
                                   sims_errors_and_bic_df100,
                                   sims_errors_and_bic_df200,
                                   sims_errors_and_bic_df500,
                                   sims_errors_and_bic_df1000)

ESTS_NAMES = c("Pearson", "FM", "Glasso", "LW","IVE", "SCE", "WSCE")#c("Pearson","LW","Sparse","FM","hatSigma0","hatSigma","WSCE")
COLVEC = c("brown","grey","pink","pink3","beige","orange2","darkorange2")
names(sims_errors_and_bic_dffull)[2:8] = ESTS_NAMES
melted_df = melt(sims_errors_and_bic_dffull)
names(melted_df)
melted_df[,3:4] = 0
for(s in 0:30){
  melted_df[(1+s*280):((s+1)*280),3:4] = melted_df[1:280,1:2]
}
melted_df = melted_df[(melted_df$variable=="Pearson")|
                        (melted_df$variable=="FM")|
                        (melted_df$variable=="Glasso")|
                        (melted_df$variable=="LW")|
                        (melted_df$variable=="IVE")|
                        (melted_df$variable=="SCE")|
                        (melted_df$variable=="WSCE"),]
melted_df$V4 = c(" ",melted_df$V4)[-1]
# reverse the order
degrees_of_freedom_labels = unique(melted_df$V4)
#reversedlabels = paste0("0",10000-as.numeric(degrees_of_freedom_labels))
reversedlabels = c(paste0("00",as.numeric(degrees_of_freedom_labels)[as.numeric(degrees_of_freedom_labels)<100]),
                   paste0("0",as.numeric(degrees_of_freedom_labels)[(as.numeric(degrees_of_freedom_labels)>=100)&(as.numeric(degrees_of_freedom_labels)<1000)]),
                   paste0(as.numeric(degrees_of_freedom_labels)[as.numeric(degrees_of_freedom_labels)>=1000]))
for(i in seq_along(degrees_of_freedom_labels)){
  melted_df$V4[melted_df$V4==degrees_of_freedom_labels[i]] = reversedlabels[i]
}
names(degrees_of_freedom_labels) = reversedlabels
plot_without_lines = ggplot(melted_df,aes(x=V4,y=value,fill=variable)) +
  geom_boxplot(color = "grey10",linewidth = .3,outlier.shape = NA) + 
  xlab("degrees of freedom") + ylab("MAE") +
  scale_fill_manual(values = c("brown","grey","pink","pink3","darkorange2","beige","orange2"),
                    labels =  expression("Pearson","FM","Glasso","LW","IVE","SCE","WSCE"),
                    name="") +
  scale_x_discrete(labels = degrees_of_freedom_labels) +
  theme(legend.position = "bottom") +
  guides(fill = guide_legend(nrow = 1))

for(i in (1:10)){
  plot_without_lines = plot_without_lines + geom_vline(xintercept = i+.5,
                                                       linetype ="dashed")
}

ggsave(plot_without_lines+theme(text = element_text(size = 16)), 
       file = "atelier/MAE_tdistribution.pdf", width = 7, height = 5)
## END model misspecification ##