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

my_theme <- theme_bw() +
  theme(strip.background = element_rect(fill = "white"), text = element_text(face="bold", size=12),
  )
theme_set(my_theme)

### Simulation 2 ###

## reading in the data ##

source("functions/get_data_covar.R")

read_plot = read_plot_FITcomps_std(filename="../Data/TFR_pieces_202311/standardized_residuals_202311/FITcomps_std_residuals_sample%i_202311.txt") # Initializing and plotting the real data 
FITcomps_std_total = read_plot$FITcomps_std_total
FITcomps_std = read_plot$FITcomps_std

read_names = read_names_FITcomps_std_total(FITcomps_std_total, covar, model="all_values")
names_by_id = read_names$names_by_id
all_min = read_names$all_min

preproc_res = preproc_FITcomps_std(all_min, names_by_id, FITcomps_std, covar)
matList_final = preproc_res$matList_final
id_min = preproc_res$id_min

## End reading the data ##

matList2 = read_matList(filename = paste(data_source,"sim_02_matList.csv",sep=""))
(sim_02_true_param = read_param(filename=paste(data_source,"sim_02_true_param.csv",sep="")))
# plot_param(matList2, filename_true_param=paste(data_source,"sim_02_true_param.csv",sep=""),
#            filename_param_fit=paste(data_source,"sim_02_param_fit.csv",sep=""), id_min=id_min)

Sigma = CovMat_03(as.matrix(sim_02_true_param), matList2, id_min=id_min)$Sigma
ests = read_ests(filename=paste(data_source,"sim_02_ests.csv",sep=""))

# known means and variances
plot_cov(matList2,Sigma,
         SigmaHat_list=read_ests(filename=paste(data_source,"sim_02_ests.csv",sep="")),
         colvec=c("brown","grey","pink","pink3","beige","orange2","darkorange2"),
         model="corY", ests_names=ESTS_NAMES,order=c(1, 2, 3, 4, 7, 5, 6))
ggsave("atelier/sim_02_ests.jpeg", width=5.3,height=4.07,device="jpeg")

# unknown means and variances
plot_cov(matList2,Sigma,
         SigmaHat_list=read_ests(filename=paste(data_source,"sim_02_ests_musigma_unknown.csv",sep="")),
         colvec=c("brown","grey","pink","pink3","beige","orange2","darkorange2"),
         model="corY", ests_names=ESTS_NAMES,order=c(1, 2, 3, 4, 7, 5, 6))
ggsave("atelier/sim_02_ests_musigmaunknown.jpeg", width=5.3,height=4.07,device="jpeg")

# heatmap of the correlation matrix
plot_heatmaps(matList2, Sigma, 
              filename=paste(data_source,"sim_02_all_matrices_plot.jpeg",sep=""))
ggsave("atelier/sim_02_all_matrices_plot.pdf", width=10.6,height=8.14)

# boxplot of error measures 
sims_errors_and_bic = read.csv(file=paste(data_source,"sim_02_sims_errors_and_bic.csv",sep=""))

param_pos = sapply(names(sims_errors_and_bic), function(s) grepl("param",s))
sims_params = sims_errors_and_bic[,param_pos]

params_1_pos = sapply(names(sims_params), function(s) grepl("param1",s))
sims_params1 = sims_params[,params_1_pos]
names(sims_params1) = PARAM1_NAMES
ggplot(melt(sims_params1),aes(x=variable,y=value)) + geom_boxplot() + 
  geom_abline(slope=0, intercept=as.numeric(sim_02_true_param[1])) +
  geom_abline(slope=0, intercept=as.numeric(sim_02_true_param[2])) +
  geom_abline(slope=0, intercept=as.numeric(sim_02_true_param[3])) +
  geom_abline(slope=0, intercept=as.numeric(sim_02_true_param[4])) +
  geom_abline(slope=0, intercept=as.numeric(sim_02_true_param[5]))

sims_errors_and_bic = sims_errors_and_bic[,!param_pos]
plot_sims(sims_errors_and_bic=sims_errors_and_bic,filename="atelier/sim_02_error_measures.pdf")

#mu and sigma unknown
# boxplot of error measures 
sims_errors_and_bic = read.csv(file=paste(data_source,"sim_02_sims_errors_and_bic_musigma_unknown.csv",sep=""))

param_pos = sapply(names(sims_errors_and_bic), function(s) grepl("param",s))
sims_params = sims_errors_and_bic[,param_pos]

params_1_pos = sapply(names(sims_params), function(s) grepl("param1",s))
sims_params1 = sims_params[,params_1_pos]
names(sims_params1) = PARAM1_NAMES
ggplot(melt(sims_params1),aes(x=variable,y=value)) + geom_boxplot() + 
  geom_abline(slope=0, intercept=as.numeric(sim_02_true_param[1])) +
  geom_abline(slope=0, intercept=as.numeric(sim_02_true_param[2])) +
  geom_abline(slope=0, intercept=as.numeric(sim_02_true_param[3])) +
  geom_abline(slope=0, intercept=as.numeric(sim_02_true_param[4])) +
  geom_abline(slope=0, intercept=as.numeric(sim_02_true_param[5]))

sims_errors_and_bic = sims_errors_and_bic[,!param_pos]
plot_sims(sims_errors_and_bic=sims_errors_and_bic,filename="atelier/sim_02_error_measures_musigma_unknown.pdf")


### End Simulation 2 ###

# varying n

my_theme <- theme_bw() +
  theme(strip.background = element_rect(fill = "white"), text = element_text(face="bold", size=12),
        axis.ticks.x=element_blank()
  )
theme_set(my_theme)

df=read.csv(file=paste(data_source,"sim_02_different_n.csv",sep=""),row.names=1)

#df1 = df[(df$est=="hatSigma")&(df$type =="MAE"),]
df$lambda=c(" ",(df$lambda))[-1]
df$lambda[(df$lambda=="14")|(df$lambda=="32")|(df$lambda=="65")] = paste("0",df$lambda[(df$lambda=="14")|(df$lambda=="32")|(df$lambda=="65")],sep="")
df$d = df$lambda
df$MAE=df$value
plot_MAE = ggplot(df[(df$est=="SCE")&(df$type=="MAE"),], aes(x=d, y=MAE)) + geom_boxplot()+scale_x_discrete(labels=c("14","32","65","115","195")) #scale_fill_hue( labels = expression("IVE","LW","WSCE","Pearson","Sparse","SCE"))

#plot=arrangeGrob(plot_MAE,plot_RMSE)
plot_MAE
ggsave(plot_MAE,filename="atelier/sim_02_different_n.pdf", width=7,height=2,device="pdf")


# model misspecification

# Version 1: No missing information
my_theme <- theme_bw() +
  theme(strip.background = element_rect(fill = "white"), text = element_text(face="bold", size=12),
        axis.ticks.x=element_blank()
  )
theme_set(my_theme)

df=read.csv(file=paste(data_source,"sim_02_modelmisspec_withWSCE.csv",sep=""),row.names=1)

df1 = df[(df$est=="hatSigma")&(df$type =="MAE"),]
df$lambda=c(" ",(df$lambda))[-1]
df[(df$est=="hatSigma"),]$est = "ZSCE"
df[(df$est=="hatSigma0"),]$est = "ZIVE"
df[(df$est=="new"),]$est = "ZZWSCE"
df[(df$est=="Sparse"),]$est = "Glasso"
df[(df$est=="Pearson"),]$est = "APearson"
df$MAE=df$value
df$estimator=df$est

copy = df[df$estimator=="BestWSCE",]
copy$estimator = "ZZWSCE"
df[df$estimator=="ZZWSCE",]=copy

df = df[df$estimator!="BestWSCE",]

df[(df$estimator=="ZZWSCE"),]$estimator = "ZAWSCE"
df[(df$estimator=="ZSCE"),]$estimator = "ZCSCE"
df[(df$estimator=="ZIVE"),]$estimator = "ZBIVE"

plot_without_lines = ggplot(df, aes(x=lambda, y=MAE, fill=estimator)) + geom_boxplot(coef = 6)
for(i in (1:10)){
  plot_without_lines = plot_without_lines + geom_vline(xintercept = i+.5)
}
plot_without_lines +
  scale_fill_manual(values = c("brown","grey","pink","pink3","beige","orange2","darkorange2"),
                    labels =  expression("Pearson","FM","Glasso","LW","IVE","SCE","WSCE"), 
                    name="") + theme(legend.position = "bottom")
ggsave("atelier/sim_02_model_misspec_partial.pdf", width=7,height=4,device="pdf")

df_1 <- df
# Version 2: Missing values
my_theme <- theme_bw() +
  theme(strip.background = element_rect(fill = "white"), text = element_text(face="bold", size=12),
        axis.ticks.x=element_blank()
  )
theme_set(my_theme)

df=read.csv(file=paste(data_source,"sim_02_modelmisspec_withWSCE_missing_values.csv",sep=""),row.names=1)

df1 = df[(df$est=="hatSigma")&(df$type =="MAE"),]
df$lambda=c(" ",(df$lambda))[-1]
df[(df$est=="hatSigma"),]$est = "ZSCE"
df[(df$est=="hatSigma0"),]$est = "ZIVE"
df[(df$est=="new"),]$est = "ZZWSCE"
df[(df$est=="Sparse"),]$est = "Glasso"
df[(df$est=="Pearson"),]$est = "APearson"
df$MAE=df$value
df$estimator=df$est

copy = df[df$estimator=="BestWSCE",]
copy$estimator = "ZZWSCE"
df[df$estimator=="ZZWSCE",]=copy

df = df[df$estimator!="BestWSCE",]

df[(df$estimator=="ZZWSCE"),]$estimator = "ZAWSCE"
df[(df$estimator=="ZSCE"),]$estimator = "ZCSCE"
df[(df$estimator=="ZIVE"),]$estimator = "ZBIVE"

plot_without_lines = ggplot(df, aes(x=lambda, y=MAE, fill=estimator)) + geom_boxplot(coef = 6)
for(i in (1:10)){
  plot_without_lines = plot_without_lines + geom_vline(xintercept = i+.5)
}
plot_without_lines +
  scale_fill_manual(values = c("brown","grey","pink","pink3","beige","orange2","darkorange2"),
                    labels =  expression("Pearson","FM","Glasso","LW","IVE","SCE","WSCE"), 
                    name="") + theme(legend.position = "bottom")
ggsave("atelier/sim_02_model_misspec_partial_missing_values.pdf", width=7,height=4,device="pdf")


### New graph with both missing and not missing 

# df=read.csv(file=paste(data_source,"sim_02_modelmisspec_withWSCE06_missing_values.csv",sep=""),row.names=1) %>%
#   mutate(est = case_when(
#   est == "BestWSCE"~"hatSigma0",
#   est == "hatSigma0"~"hatSigma",
#   est == "hatSigma"~"BestWSCE",
#   TRUE~est
# ))
#read.csv(file=paste(data_source,"sim_02_modelmisspec06.csv",sep=""),row.names=1)

df_1$full ="Full"
df$full ="Missing"
df_full=rbind(df,df_1)

plot_without_lines <- df_full %>% 
   filter(est!="new"  & type =="MAE") %>% 
   ggplot( aes(x=lambda, y=value, fill=estimator)) + geom_boxplot(coef = 6) +
   facet_grid(full~.)+
   scale_fill_manual(values = c("brown","grey","pink","pink3","darkorange2","beige","orange2"),
                     labels =  expression("Pearson","FM","Glasso","LW","IVE","SCE","WSCE"),
                     name="") + theme(legend.position = "bottom") + 
   labs(y = "MAE", x = expression(lambda)) + guides(fill = guide_legend(nrow = 1))
   
 for(i in (1:10)){
   plot_without_lines = plot_without_lines + geom_vline(xintercept = i+.5, 
                                                        linetype ="dashed")
 }
 
ggsave(plot_without_lines, file = "atelier/MAE_lambda.pdf", width = 7, height = 5)
