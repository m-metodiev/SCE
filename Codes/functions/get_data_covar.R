library(tidyverse)
library(Matrix)

covar_pays <- read.table(file="../Data/covariable/DistCorPredictors.txt",sep=",", header= TRUE) %>%
  mutate(id_col = as.numeric(as.factor(iso_o)),
         id_row = as.numeric(as.factor(iso_d)))

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

  table(covar_pays_full %>% filter(comcol ==1) %>% pull(name_d))
  table(covar_pays_full %>% filter(comcol ==1) %>% pull(name_o))
  tb1 <-table(covar_pays_full$name_d, covar_pays_full$comcol)
  tb2 <-table(covar_pays_full$name_o, covar_pays_full$comcol)
  cbind(tb2[tb1[,2]!=tb2[,2],],tb1[tb1[,2]!=tb2[,2],])
  
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

tb1 <-table(covar$name_d, covar$comcol)
tb2 <-table(covar$name_o, covar$comcol)
cbind(tb2[tb1[,2]!=tb2[,2],],tb1[tb1[,2]!=tb2[,2],])

n_pays <- length(unique(covar_pays$iso_d))
Gb <- sparseMatrix(i = covar_pays$id_row, j = covar$id_col, x = covar$comcol)
Gc <- sparseMatrix(i = covar_pays$id_row, j = covar$id_col, x = covar$same_reg)
Gd <- matrix(1, ncol = n_pays,nrow = n_pays)
Ga_inv <- sparseMatrix(i = covar_pays$id_row, j = covar$id_col, x = covar$contig)
Id <- Diagonal(n_pays)