library(tidyverse)
library(Matrix)

covar_pays <- read.table(file="../Data/covariable/DistCorPredictors.txt",sep=",", header= TRUE) %>%
  mutate(id_col = as.numeric(as.factor(iso_o)),
         id_row = as.numeric(as.factor(iso_d)))


# reset the comcols
# library(cepiigeodist)
# data(dist_cepii)
# 
# test = dist_cepii
# test$numc_o = test$iso_o
# test$numc_d = test$iso_d
# 
# for(str in unique(test$iso_o)){
#   test$numc_o[dist_cepii$iso_o == str] = geo_cepii$cnum[geo_cepii$iso3 == str]
#   test$numc_d[dist_cepii$iso_d == str] = geo_cepii$cnum[geo_cepii$iso3 == str]
# }]
# test$unique_id = paste(test$numc_o,test$numc_d)
# covar_pays$unique_id = paste(covar_pays$iso_o, covar_pays$iso_d)
# test_final = left_join(covar_pays ,test, by =c("unique_id" = "unique_id"))

# dens <- function(ga, gb ,gc,gd,st){
#   S_inv <- Ga_inv  %*% solve(ga * Id + st * Ga_inv +
#                                gd * Ga_inv %*% Gd +
#                                gb * Ga_inv %*% Gb +
#                                gc * Ga_inv %*% Gc )
#   det(S_inv)^0.5 *exp(- 0.5 * t(x) %*% S_inv %*% x)
# }
# 
# 
# ## to which countries correspond the standardized residuals?
# 
# 
library(readr)

WPP2008_F01_LOCATIONS <- read_csv("../Data/TFR_pieces_202311/WPP2008_F01_LOCATIONS.csv")

countries_num <- WPP2008_F01_LOCATIONS %>%
  filter(!is.na(`ISO3 Alpha-code`))
dim(countries_num)


table(countries_num$reg_name)
# 
# ## add country name and region in covar_pays fr iso_d and iso_o
# 
# 
countries_d = countries_num %>%
  transmute(iso_d = country_code,
            reg_d = reg_code,
            name_d = country_name,
            code_d=`ISO3 Alpha-code`)
countries_o = countries_num %>%
  transmute(iso_o = country_code,
            reg_o = reg_code,
            name_o = country_name,
            code_o=`ISO3 Alpha-code`)



covar_pays_full <- covar_pays %>%
  left_join(countries_d, by = "iso_d")%>%
  left_join(countries_o, by = "iso_o") %>%
  mutate(same_reg = as.numeric(reg_d==reg_o))
countries_num %>% select(country_code, reg_name,reg_code) %>%
  spread(key=country_code, value = reg_name)
# 
# 
# ## graphe des pays adjacents 
# library(ggraph)
# library(igraph)
# 
# # Let's look at the graph of the continuous countries to check 
# edges <- covar_pays_full %>% filter(contig ==1 & ! name_o==name_d) %>% 
#   transmute(from=name_o, to = name_d)
# vertices_countries = countries_d %>% mutate(name= name_d)
# flareGraph <- graph_from_data_frame(edges)
# ggraph(flareGraph, 'igraph', algorithm = 'kk') + 
#   geom_edge_link() +
#   ggforce::theme_no_axes()+
#   geom_node_text(aes(label = name), repel=TRUE)
# 
# 
# edges <- covar_pays_full %>% filter(contig ==1 & ! name_o==name_d) %>% 
#   transmute(from=name_o, to = name_d)
# vertices_countries = countries_d %>%
#   filter(name_d %in% unique(edges$from)) %>% 
#   mutate(name= name_d) %>% 
#   select(name, reg_d)
# flareGraph <- graph_from_data_frame(edges, vertices = vertices_countries) 
# ggraph(flareGraph, 'igraph', algorithm = 'kk') + 
#   geom_edge_link() +
#   ggforce::theme_no_axes()+
#   geom_node_text(aes(label = name, color = reg_d), repel=TRUE)
# 
  table(covar_pays_full %>% filter(comcol ==1) %>% pull(name_d))
  table(covar_pays_full %>% filter(comcol ==1) %>% pull(name_o))
# 
#  
#  p_12 <- covar_pays_full %>% filter(iso_d==12 &comcol==1) %>% pull(name_o)
#  
#  p_21 <- covar_pays_full %>% filter(iso_o==12 & comcol==1)%>% pull(name_d)
#  p_12[!(p_12%in% p_21)]
# 
#  ## on a un problème avec mayotte qui ne semble pas être lié aux autres.. 
#  
#  covar_pays_full %>% filter(name_d =="Mayotte") %>% View()
  tb1 <-table(covar_pays_full$name_d, covar_pays_full$comcol)
  tb2 <-table(covar_pays_full$name_o, covar_pays_full$comcol)
  cbind(tb2[tb1[,2]!=tb2[,2],],tb1[tb1[,2]!=tb2[,2],])
  
  
  # pb Mayotte Channel Islands, United States Virgin Islands,
  # qui sont mises avec les autres quand elles sont d mais pas o.. 
  
  #Vanuatu?? étrange à voir après 
 pb <- c("Mayotte", "Channel Islands", "United States Virgin Islands")
 samecol_fun<-  function(x){covar_pays_full %>% 
   filter(name_o ==x) %>% 
   filter(comcol==1) %>% 
   pull(name_d)
 }
 samecols <- lapply(pb, samecol_fun)
 names(samecols) <- pb
 
covar <-  covar_pays_full  %>% 
   mutate(same_may =  case_when(
     name_d == "Mayotte" & name_o %in% samecols$Mayotte ~1,
     TRUE ~ 0
   ), 
   same_chan =  case_when(
     name_d == "Channel Islands" & 
       name_o %in% samecols[["Channel Islands"]] ~1,
     TRUE ~ 0
   ), 
   same_virg =  case_when(
     name_d == "United States Virgin Islands" & 
       name_o %in% samecols[["United States Virgin Islands"]] ~1,
     TRUE ~ 0
   ),
   comcol = comcol  + same_may+same_virg+same_chan)
table(covar$comcol)
sum(covar_pays_full$comcol) + length(unlist(samecols))
# ok 

## vanuatu ? ah non en faite c'est qu'il avait deux colons 

tb1 <-table(covar$name_d, covar$comcol)
tb2 <-table(covar$name_o, covar$comcol)
cbind(tb2[tb1[,2]!=tb2[,2],],tb1[tb1[,2]!=tb2[,2],])

# tout bon ! 

n_pays <- length(unique(covar_pays$iso_d))
Gb <- sparseMatrix(i = covar_pays$id_row, j = covar$id_col, x = covar$comcol)
Gc <- sparseMatrix(i = covar_pays$id_row, j = covar$id_col, x = covar$same_reg)
Gd <- matrix(1, ncol = n_pays,nrow = n_pays)
Ga_inv <- sparseMatrix(i = covar_pays$id_row, j = covar$id_col, x = covar$contig)
Id <- Diagonal(n_pays)


