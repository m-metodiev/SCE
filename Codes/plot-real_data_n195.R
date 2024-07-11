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
ESTS_NAMES = c("Pearson","LW","Sparse","hatSigma0","hatSigma")

### True data 2: n=195 ###

## read data ##

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

## End read data ##

# #Comparing to Adrian's results:
# plot_param(matList_final, filename_true_param="sim_02_true_param.csv",
#            filename_param_fit="sim_final_n195_param_fit.csv", id_min=id_min,
#            error=FALSE)
# plot_param(matList_final, filename_true_param="sim_02_true_param.csv",
#            filename_param_fit="sim_final_n195_combined_param_fit.csv", id_min=id_min,
#            error=FALSE)
# 
# plot_param(matList_final, filename_true_param=NULL,
#            filename_param_fit="sim_final_n195_param_fit.csv", id_min=id_min)
# 
# plot_param(matList_final, filename_true_param=NULL,
#            filename_param_fit="sim_final_n74_combined_param_fit.csv", id_min=id_min)

#Comparing to Adrian's results:
# plot_param(matList_final, filename_true_param=paste(data_source,"sim_02_true_param.csv",sep=""),
#            filename_param_fit=paste(data_source,"sim_final_n195_param_fit.csv",sep=""), 
#            id_min=id_min,
#            error=FALSE)
# plot_param(matList_final, filename_true_param=paste(data_source,"sim_02_true_param.csv",sep=""),
#            filename_param_fit=paste(data_source,"sim_final_n195_combined_param_fit.csv",sep=""),
#            id_min=id_min,
#            error=FALSE)
# plot_param(matList_final, filename_true_param=paste(data_source,"sim_02_true_param.csv",sep=""),
#            filename_param_fit=paste(data_source,"sim_final_n195_param_fit.csv",sep=""), 
#            id_min=id_min,
#            error=FALSE,version="scatterplot")
# plot_param(matList_final, filename_true_param=paste(data_source,"sim_02_true_param.csv",sep=""),
#            filename_param_fit=paste(data_source,"sim_final_n195_combined_param_fit.csv",sep=""), 
#            id_min=id_min,
#            error=FALSE,version="scatterplot")


## Compare results
Sigma_FosdickRaftery = as.matrix(CovMat_03(parm=c(.05,.09,.11,.26),matList = matList_final,
                                           id_min=id_min,combined_effects = "FosdickRaftery")$Sigma)
Sigma_base_model = as.matrix(CovMat_03(parm=as.matrix(read_param(filename=paste(data_source,"sim_final_n195_param_fit.csv",sep=""))),
                                       matList=matList_final,
                                       id_min=id_min)$Sigma)
Sigma_interation_model = as.matrix(CovMat_03(parm=as.matrix(read_param(filename=paste(data_source,"sim_final_n195_combined_param_fit.csv",sep=""))),
                                             matList=matList_final,
                                             id_min=id_min, link=combined_matList)$Sigma)
# plot_heatmaps(matList_final, Sigma=Sigma_FosdickRaftery, 
#               Sigma2=Sigma_base_model, Sigma3=Sigma_interation_model, 
#               filename=paste(data_source,"sim_final_n195_all_matrices_plot.jpeg",sep=""))

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
  mutate(short_colonizer1 = map_chr(iso3, ~get_count_info(., missing_iso, covar, geo_cepii))) %>% 
  bind_rows(geo_cepii %>% select(iso3,short_colonizer1))

## Now let's merge everything we need 

info_coutries <- countries_num %>% 
  select(country_name,`ISO3 Alpha-code`, reg_name) %>%
  left_join(colon_full, by =c("ISO3 Alpha-code" = "iso3"))  %>% 
  unique()

#Y = FITcomps_std_total[2:12,which(all_min==1)]
# Y = Y[,sapply(id_min,function(id) which(preproc_res$FITcomps_std_iso==preproc_res$iso_id_key[id]))]
# n=dim(Y)[2]
# imputed = imputePCA(Y,ncp=4)$completeObs
# LWest = linearShrinkLWEst(imputed)#res_list[[i]][[3]][[2]]
SigmaHatSCE = Sigma_interation_model
# cov_est = cov(imputed)
# testfunc = function(lamb) mean(abs(((1-lamb)*diag(n)+lamb*(cov_est))-LWest))
# lambda=which.min(sapply((1:10000)/10000,testfunc))/10000

ests = read_ests(filename=paste(data_source,"sim_final_n195_combined_ests.csv",sep=""))

main_list = plot_heatmaps(matList_final, 
              Sigma=ests[[7]],
              filename="atelier/sim_final_n195_full_data_covmatWSCE.jpeg",
              show_regions=TRUE, Sigma_true=SigmaHatSCE,
              names_by_id=names_by_id, id_min=id_min,
              iso3=iso3, info_coutries=info_coutries)

# cl3 = plot_heatmaps(matList_final, 
#                     Sigma=ests[[1]],
#                     filename="atelier/sim_final_n195_full_data_covmatPEARSON.jpeg",
#                     show_regions=TRUE, Sigma_true=SigmaHatSCE,
#                     names_by_id=names_by_id, id_min=id_min,
#                     iso3=iso3, info_coutries=info_coutries)

plot_heatmaps(matList_final, 
              Sigma=Sigma_base_model,
              filename="atelier/sim_final_n195_full_data_covmatbasemodel.jpeg",
              show_regions=TRUE, Sigma_true=SigmaHatSCE,
              names_by_id=names_by_id, id_min=id_min,
              iso3=iso3, info_coutries=info_coutries)

matList_final$Al[[1]][is.na(matList_final$Al[[1]])]=0


cl2 = plot_heatmaps(matList_final, 
              Sigma=SigmaHatSCE,
              filename="atelier/sim_final_n195_full_data_covmat.jpeg",
              show_regions=TRUE, Sigma_cluster=ests[[7]],
              names_by_id=names_by_id, id_min=id_min,
              iso3=iso3, info_coutries=info_coutries,
              main_range = main_list$main_range, main_color = main_list$main_color)

#TODO: remove this
cl2 = plot_heatmaps(matList_final, 
                    Sigma=SigmaHatSCE,
                    filename="atelier/sim_final_n195_full_data_covmat.jpeg",
                    show_regions=TRUE, Sigma_true = Sigma_base_model,
                    names_by_id=names_by_id, id_min=id_min,
                    iso3=iso3, info_coutries=info_coutries)

adjustedRandIndex(cl1,cl2)
adjustedRandIndex(cl1,cl3)

### End True data 2 ###

