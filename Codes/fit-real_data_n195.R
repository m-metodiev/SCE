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

### START: check identifiability (this may take a long time) ###

# grid-search for rho/=rho' 
# (independence of the matrices is checked automatically)
beta=1.5
num_grid_sim=100
xi = (1:(num_grid_sim+1))/(num_grid_sim+1)
#tan-hyperbolic-spaced grid because rho approaches 1
rho_vec = (1-tanh(beta*(1+xi))/tanh(beta))/(min((1-tanh(beta*(1+xi))/tanh(beta))))
rho_vec = rho_vec[-length(rho_vec)]

library(lpSolve)

for(rho_1 in rho_vec){
  for(rho_2 in rho_vec){
    if(rho_1!=rho_2){
      K = length(matList_final$Fk)
      f.obj = c(rep(0,K),1,rep(0,K),1)
      
      
      
      G_matrix_1 = calc_tilde_G_inv(M=matList_final$Ml[[1]],A=matList_final$Al[[1]],
                                    rho=rho_1)[id_min,id_min]
      G_matrix_2 = calc_tilde_G_inv(M=matList_final$Ml[[1]],A=matList_final$Al[[1]],
                                    rho=rho_2)[id_min,id_min]
      
      for(k in 1:K){
        diag(matList_final$Fk[[k]])=1
      }
      
      # matrix combination constraints
      f.con = cbind(sapply(1:K,function(k) c(as.matrix(matList_final$Fk[[k]]))),
                    c(as.matrix(G_matrix_1)),
                    -sapply(1:K,function(k) c(as.matrix(matList_final$Fk[[k]]))),
                    -c(as.matrix(G_matrix_2)))
      f.rhs = rep(0,dim(f.con)[1])
      
      # simplex constraints
      f.con = rbind(f.con, c(rep(1,K+1),rep(0,K+1)))
      f.con = rbind(f.con, c(rep(0,K+1),rep(1,K+1)))
      f.dir = rep("=", dim(f.con)[1])
      f.rhs = c(f.rhs,1,1)
      
      # positivity constraints
      f.con = rbind(f.con, diag(2*K+2))
      f.dir = c(f.dir,rep(">=", 2*K+2))
      f.rhs = c(f.rhs, rep(0,2*K+2))
      
      solution <- lp("max",f.obj,f.con,f.dir,f.rhs)
      if(solution$objval!=0){
        print("WARNING: NOT IDENTIFIABLE")
        print(solution$solution)
      }
    }
  }
}

### END: check identifiability ###

fit_param(n, p, matList_final, sim_final, id_min=id_min,
          filename_param_fit=paste(data_source,"sim_final_n195_param_fit.csv",sep=""),
          filename_ests = paste(data_source,"sim_final_n195_ests.csv",sep=""),
          filename_bic = paste(data_source,"sim_final_n195_bic.csv",sep=""))
read.csv(paste(data_source,"sim_final_n195_param_fit.csv",sep=""))

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
  
  combined_matList_partial_der = function(matList, link_matList, tilde_G_inv_partial_rho){
    #browser()
    link_matList=c(matList$Fk,matList$Gl)
    n = dim(link_matList[[1]])[1]
    matList_full = c(matList$Fk,matList$Gl)
    counter = length(matList_full)
    sequence = seq_along(matList_full)
    #calculate all possible Hadamard-products; Exclude global effect matrix
    for(i in sequence){
      matList_full[[i]] = matrix(0,n,n) # does not depend on rho
    }
    
    matList_full[[length(matList$Fk)+1]] = tilde_G_inv_partial_rho
    #browser()
    for(i in sequence[-length(matList$Fk)]){
      for(j in sequence[c(-(1:i),-length(matList$Fk))]){
        
        counter = counter + 1
        
        # use product rule whenever tilde_G_inv is included
        if(i==(length(matList$Fk)+1)){
          matList_full[[counter]] = tilde_G_inv_partial_rho * link_matList[[j]]
        } else if(j==(length(matList$Fk)+1)){
          matList_full[[counter]] = link_matList[[i]] * tilde_G_inv_partial_rho
        } else{
          matList_full[[counter]] = matrix(0,n,n)
        }
      }
    }
    #browser()
    # Only select relevant matrices
    comb_mat = matList_full
    mat_list_full = list()
    counter=1
    for(j in which(as.numeric(intToBits(k))[1:7]==1)){
      mat_list_full[[counter]] = comb_mat[[j]]
      counter = counter+1
    }
    return(mat_list_full)
  }
  #browser()
  this0 = fit_param(n, p, matList_final, sim_final, id_min=id_min,
                    link=combined_matList_partial,
                    link_der_rho = combined_matList_partial_der, save=FALSE)
  #browser()
  this_list[[k]] = this0
  print(k)
  #browser()
}

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

# When the spatial effect is missing, we are down TWO parameters, not 1
p=11
sort_bics = sort(bics)[simple_sequence]
sort_bics[is.na(df$'3contig')] = sort_bics[is.na(df$'3contig')] - log(p)
rownames(df) = round(sort_bics,1)
write.csv(df,file=paste(data_source,"sim_final_n195_model_choice.csv",sep=""))

round(df[1,],6)
matList_final$Gl[[1]] = matrix(0,n,n)
write_param(t(c(df_test[1,])),filename=paste(data_source,"sim_final_n195_combined_param_fit.csv",sep=""),
            matList = matList_final, link=combined_matList)
round(df[rownames(df)=="0",],3)

# calculate SCE on final model
combined_matList_partial = function(matList){
  #browser()
  comb_mat = combined_matList(matList)
  mat_list_full = list()
  counter=1
  for(j in which(c(1,1,1,1,0,1,1)==1)){
    mat_list_full[[counter]] = comb_mat[[j]]
    counter = counter+1
  }
  return(mat_list_full)
}

combined_matList_partial_der = function(matList, link_matList, tilde_G_inv_partial_rho){
  #browser()
  link_matList=c(matList$Fk,matList$Gl)
  n = dim(link_matList[[1]])[1]
  matList_full = c(matList$Fk,matList$Gl)
  counter = length(matList_full)
  sequence = seq_along(matList_full)
  #calculate all possible Hadamard-products; Exclude global effect matrix
  for(i in sequence){
    matList_full[[i]] = matrix(0,n,n) # does not depend on rho
  }
  
  matList_full[[length(matList$Fk)+1]] = tilde_G_inv_partial_rho
  #browser()
  for(i in sequence[-length(matList$Fk)]){
    for(j in sequence[c(-(1:i),-length(matList$Fk))]){
      
      counter = counter + 1
      
      # use product rule whenever tilde_G_inv is included
      if(i==(length(matList$Fk)+1)){
        matList_full[[counter]] = tilde_G_inv_partial_rho * link_matList[[j]]
      } else if(j==(length(matList$Fk)+1)){
        matList_full[[counter]] = link_matList[[i]] * tilde_G_inv_partial_rho
      } else{
        matList_full[[counter]] = matrix(0,n,n)
      }
    }
  }
  #browser()
  # Only select relevant matrices
  comb_mat = matList_full
  mat_list_full = list()
  counter=1
  for(j in which(c(1,1,1,1,0,1,1)==1)){
    mat_list_full[[counter]] = comb_mat[[j]]
    counter = counter+1
  }
  return(mat_list_full)
}

fit_param(n, p, matList_final, sim_final, id_min=id_min,
          filename_param_fit=paste(data_source,"sim_final_n195_combined_param_fit.csv",sep=""),
          filename_ests = paste(data_source,"sim_final_n195_combined_ests.csv",sep=""),
          filename_bic = paste(data_source,"sim_final_n195_combined_bic.csv",sep=""),
          link=combined_matList_partial, 
          link_der_rho = combined_matList_partial_der, compute_WSCE=TRUE)
read.csv(paste(data_source,"sim_final_n195_combined_param_fit.csv",sep=""))
# the WSCE was equal to the SCE
fit_param(n, p, matList_final, sim_final, id_min=id_min,
          filename_param_fit=paste(data_source,"sim_final_n195_combined_param_fit.csv",sep=""),
          filename_ests = paste(data_source,"sim_final_n195_combined_ests.csv",sep=""),
          filename_bic = paste(data_source,"sim_final_n195_combined_bic.csv",sep=""),
          link=combined_matList_partial, link_der_rho = combined_matList_partial_der)
parm=read.csv(paste(data_source,"sim_final_n195_combined_param_fit.csv",sep=""))

covMatstuff = CovMat_03(parm=as.numeric(c(as.matrix(parm))[-1]) ,matList=matList_final,id_min=id_min,link=combined_matList)
ml_combined = covMatstuff$matList_combined
alpha_beta = covMatstuff$alpha_beta

sapply(seq_along(alpha_beta), function(i) alpha_beta[i]*sum(ml_combined_supp[[i]]*ml_combined[[i]])/sum(ml_combined_supp[[i]]))

parm=read.csv(paste(data_source,"sim_final_n195_param_fit.csv",sep=""))

covMatstuff = CovMat_03(parm=c(as.numeric(c(as.matrix(parm))[-1])[-5],0,0,0,as.numeric(c(as.matrix(parm))[-1])[5]) ,matList=matList_final,id_min=id_min,link=combined_matList)
ml_combined = covMatstuff$matList_combined
alpha_beta = covMatstuff$alpha_beta

sapply(seq_along(alpha_beta), function(i) alpha_beta[i]*sum(ml_combined_supp[[i]]*ml_combined[[i]])/sum(ml_combined_supp[[i]]))
