rm(list=ls())
source("functions/FunctionsCovTFR_02.R")
source("functions/funcs_frob_norm_opt.R")
source("functions/cov_TFR_sim_funcs.R")
source("functions/cov_TFR_fit_funcs.R")
source("functions/cov_TFR_plot_funcs.R")
source("functions/cov_TFR_data_funcs.R")
data_source = "data/"

## Path 2: Use ALL countries

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

seed=3
set.seed(seed)

Y = FITcomps_std_total[2:12,which(all_min==1)]
Y = Y[,sapply(id_min,function(id) which(preproc_res$FITcomps_std_iso==preproc_res$iso_id_key[id]))]
n=dim(Y)[2]
corY = matrix(ncol=n,nrow=n)

#use pairwise correlation estimates
for(i in (1:n)){
  for(j in (i:n)){
    #browser()
    Yi_notmissing = !is.na(Y[,i])
    Yj_notmissing = !is.na(Y[,j])
    corY_ij = cor(Y[Yi_notmissing&Yj_notmissing,i], Y[Yi_notmissing&Yj_notmissing,j])
    corY[i,j] = corY_ij
    corY[j,i] = corY_ij
  }
}
mean(is.na(corY))
#library(missMDA)
#Y_nonmissing = imputePCA(Y,ncp=6)$completeObs
#corY = cor(Y_nonmissing)
n = dim(matList_final$Fk[[1]])[1]
p=11
sim_final = list(corY=corY,Y=Y)

source("functions/cov_TFR_fit_funcs.R")
fit_param(n, p, matList_final, sim_final, id_min=id_min,
          filename_param_fit=paste(data_source,"sim_final_n195_param_fit.csv",sep=""),
          filename_ests = paste(data_source,"sim_final_n195_ests.csv",sep=""),
          filename_bic = paste(data_source,"sim_final_n195_bic.csv",sep=""))
read.csv(paste(data_source,"sim_final_n195_param_fit.csv",sep=""))

#TODO: Exclude 1 interaction effect

fit_param(n, p, matList_final, sim_final, id_min=id_min,
          filename_param_fit=paste(data_source,"sim_final_n195_combined_param_fit.csv",sep=""),
          filename_ests = paste(data_source,"sim_final_n195_combined_ests.csv",sep=""),
          filename_bic = paste(data_source,"sim_final_n195_combined_bic.csv",sep=""),
          link=combined_matList, compute_WSCE=TRUE)
read.csv(paste(data_source,"sim_final_n195_combined_param_fit.csv",sep=""))

this_list = list()
this = c(0,0,0,0,0,0,0)
for(k in (1:(2^length(this)))){
  combined_matList_partial = function(matList){
    #browser()
    comb_mat = combined_matList(matList)
    mat_list_full = list()
    counter=1
    for(j in which(as.numeric(intToBits(k))[1:7]==1)){
      mat_list_full[[counter]] = comb_mat[[j]]
      counter = counter+1
    }
    return(mat_list_full)
  }
  this0= fit_param(n, p, matList_final, sim_final, id_min=id_min,
                  link=combined_matList_partial, has_missingvalues = TRUE,
                  save=FALSE)
  #browser()
  this_list[[k]] = this0
  #browser()
}

save(this_list,file="huge_data_10.Rdata")
load(file="huge_data_10.Rdata")
bics = sapply(1:length(this_list),function(s) this_list[[s]][[3]])
simple_model_bic = sort(bics)[(sort(bics)<4796.5)&(sort(bics)>4796)]
simple_model_ix = which((sort(bics)<4796.5)&(sort(bics)>4796))
bics = bics - simple_model_bic

plot(sort(bics),xlab="Index of the model",ylab="BIC",main="sorted BIC values")
ix_optim_model = which.min(bics)
model_params=as.numeric(intToBits(ix_optim_model))[1:7]
model_params[model_params==1]=this_list[[ix_optim_model]][[1]][-length(this_list[[ix_optim_model]])]
model_params = c(model_params,this_list[[ix_optim_model]][[1]][length(this_list[[ix_optim_model]][[1]])])


model_param_matrix = matrix(nrow=length(this_list),ncol=length(model_params))
counter=1
for(i in sort(bics,index.return=TRUE)$ix){
  params = as.numeric(intToBits(i))[1:7]
  params[params==1] = this_list[[i]][[1]][-length(this_list[[i]][[1]])]
  params = c(params,this_list[[i]][[1]][length(this_list[[i]][[1]])])
  #browser()
  model_param_matrix[counter,] = params
  counter=counter+1
}

# sequence which does not include cases where the regional effect is excluded,
# but the combined effect is included (same for the comcol effect)
simple_sequence = which((((((model_param_matrix[,1]==0) & (model_param_matrix[,6]!=0))|
                            (((model_param_matrix[,2]==0) & (model_param_matrix[,7]!=0)))|
                            (((model_param_matrix[,1]==0) & (model_param_matrix[,5]!=0)))|
                             (((model_param_matrix[,2]==0) & (model_param_matrix[,5]!=0)))|
                             (((model_param_matrix[,4]==0) & (model_param_matrix[,6]!=0)))|
                             (((model_param_matrix[,4]==0) & (model_param_matrix[,7]!=0))))+1)%%2)==1)
length(simple_sequence)

#Determine support of the matrices
matList_supp = matList_final

#Set to 1 if countries are not neighbors
matList_supp$Al[[1]][is.na(matList_supp$Al[[1]])] = 0
matList_supp$Gl = (matList_supp$Al[[1]][id_min,id_min] != 0) + 0

ml_combined_supp = combined_matList(matList_supp) 
# the support of the product is the product of the supports

#Calculate average effect (=mean effect over matrix support)
avg_effects = function(parm){
  covMatstuff = CovMat_03(parm=parm ,matList=matList_final,id_min=id_min,link=combined_matList)
  ml_combined = covMatstuff$matList_combined
  alpha_beta = covMatstuff$alpha_beta
  
  sapply(seq_along(alpha_beta), function(i) alpha_beta[i]*sum(ml_combined_supp[[i]]*ml_combined[[i]])/sum(ml_combined_supp[[i]]))
}

df_test = model_param_matrix[simple_sequence,]
df = as.data.frame(t(apply(df_test,1,avg_effects)))
df[df==0] = NA
names(df) = c("0comcol","1reg", "2global", "3contig", 
              "4comcol and reg", "5comcol and contig", 
              "6reg and contig")

# correct: When the spatial effect is missing, we are down TWO parameters, not 1
p=11
sort_bics = sort(bics)[simple_sequence]
sort_bics[is.na(df$'3contig')] = sort_bics[is.na(df$'3contig')] - log(p)
rownames(df) = round(sort_bics,1)
library(viridis)
ggheatmap(df) +
  scale_fill_viridis() +
  guides(fill=guide_colorbar("Value")) +
  scale_x_discrete( labels = expression("comcol","sameRegion","intercept","contig","comcol and sameRegion","comcol and contig", "sameRegion and contig"))+
  scale_y_discrete(limits = rownames(df)) +
  annotate('rect', xmin = 0.5, xmax = 7.5, ymin = 0.5, ymax = 1.5,
           fill = NA, color = 'magenta', size = 1)
round(df[1,],5)
matList_final$Gl[[1]] = matrix(0,n,n)
write_param(t(c(df_test[1,])),filename=paste(data_source,"sim_final_n195_combined_param_fit.csv",sep=""),
            matList = matList_final, link=combined_matList)
ggsave("atelier/real_data_n195_modelchoice.pdf", width=12.6,height=7.14)
round(df[rownames(df)=="0",],3)


combined_matList = function(matList){
  #browser()
  matList_full = c(matList$Fk,matList$Gl)
  counter = length(matList_full)
  sequence = seq_along(matList_full)
  #calculate all possible Hadamard-products; Exclude global effect matrix
  for(i in sequence[-length(matList$Fk)]){
    for(j in sequence[c(-(1:i),-length(matList$Fk))]){
      counter = counter + 1
      matList_full[[counter]] = matList_full[[i]] * matList_full[[j]]
    }
  }
  #browser()
  return(matList_full)
}

read.csv(paste(data_source,"sim_final_n195_combined_bic.csv",sep=""))$bic-
  read.csv(paste(data_source,"sim_final_n195_bic.csv",sep=""))$bic
