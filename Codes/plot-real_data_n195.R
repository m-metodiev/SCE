rm(list=ls())
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
ESTS_NAMES = c("Pearson","LW","Sparse","hatSigma0","hatSigma")

df=read.csv(file=paste(data_source,"sim_final_n195_model_choice.csv",sep=""))
rownames(df)=df$X
df$X=NULL
library(viridis)
ggheatmap(df) +
  scale_fill_viridis(na.value = "lightgrey") +
  guides(fill=guide_colorbar("Value")) +
  scale_x_discrete( labels = expression("comcol","sameRegion","intercept","contig","comcol & sameRegion","comcol & contig", "sameRegion & contig"))+
  scale_y_discrete(limits = rownames(df)) +
  annotate('rect', xmin = 0.5, xmax = 7.5, ymin = 0.5, ymax = 1.5,
           fill = NA, color = 'magenta', size = 1) +
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=16))
ggsave("atelier/real_data_n195_modelchoice.pdf", width=16,height=10)

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

Sigma_base_model = as.matrix(CovMat_03(parm=as.matrix(read_param(filename=paste(data_source,"sim_final_n195_param_fit.csv",sep=""))),
                                       matList=matList_final,
                                       id_min=id_min)$Sigma)

Sigma_interation_model = as.matrix(CovMat_03(parm=as.matrix(read_param(filename=paste(data_source,"sim_final_n195_combined_param_fit.csv",sep=""))),
                                             matList=matList_final,
                                             id_min=id_min, link=combined_matList)$Sigma)

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

SigmaHatSCE = Sigma_interation_model

ests = read_ests(filename=paste(data_source,"sim_final_n195_combined_ests.csv",sep=""))


plot_heatmaps(matList_final, 
              Sigma=Sigma_base_model,
              filename="atelier/sim_final_n195_full_data_covmatbasemodel.jpeg",
              show_regions=TRUE, Sigma_true=SigmaHatSCE,Sigma_cluster=Sigma_base_model,
              names_by_id=names_by_id, id_min=id_min,
              iso3=iso3, info_coutries=info_coutries)

matList_final$Al[[1]][is.na(matList_final$Al[[1]])]=0


plot_heatmaps(matList_final, 
                    Sigma=SigmaHatSCE,
                    filename="atelier/sim_final_n195_full_data_covmat.jpeg",
                    show_regions=TRUE,
                    names_by_id=names_by_id, id_min=id_min,
                    iso3=iso3, info_coutries=info_coutries)
