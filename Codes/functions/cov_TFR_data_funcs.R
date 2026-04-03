write_param = function(param_fit1,filename, matList, interaction_effects=list()){
  #browser()
  param_fit_csv = as.data.frame(t(c(param_fit1)))
  rownames(param_fit_csv) = "param"
  if(length(combined_matList(matList,interaction_effects=interaction_effects)$matList_full)==4){
    colnames(param_fit_csv) = c("comcol","reg", "global", "contig beta", "contig beta")
  } else{
    colnames(param_fit_csv) = c("comcol","reg", "global", "contig beta",
                                "comcol and reg", "comcol and config",
                                "reg and config" ,"contig beta")
  }
  write.csv(param_fit_csv, file=filename)
}
read_param = function(filename){
  return(read.csv(filename,row.names=1))
}
write_matList = function(matList,filename){
  matList$Al[[1]][is.na(matList$Al[[1]])]=0
  
  #set values "outside of the matrix" to NA
  s = dim(matList$Al[[1]])[1]
  res = matrix(nrow=s*s,ncol=length(c(matList$Ml,matList$Al,matList$Fk)))
  for(i in seq_along(res[1,])){
    s = dim(matList$Al[[1]])[1]
    mat_i = as.matrix(c(matList$Ml,matList$Al,matList$Fk)[[i]])
    dimension = dim(mat_i)[1]
    res[,i][1:(dimension*dimension)] = c(mat_i)
  }
  write.csv(as.data.frame(res), file=filename)
}
read_matList = function(filename){
  full_matrix = as.matrix(read.csv(file=filename))
  full_matrix = full_matrix[,2:dim(full_matrix)[2]]
  s = round(sqrt(length(full_matrix[,1])))
  matList = list(Fk=list(),Gl=list(),
                 Ml=list(matrix(full_matrix[,1],s,s)),
                 Al=list(matrix(full_matrix[,2],s,s)))
  for(i in 3:dim(full_matrix)[2]){
    dimension = round(sqrt(sum(!is.na(full_matrix[,i]))))
    matList$Fk[[i-2]] = matrix(full_matrix[1:(dimension*dimension),i],dimension,dimension)
  }
  return(matList)
}
read_ests = function(filename){
  full_matrix = as.matrix(read.csv(file=filename))
  full_matrix = full_matrix[,2:dim(full_matrix)[2]]
  dimension = round(sqrt(length(full_matrix[,1])))
  ests = list()
  for(i in (1:dim(full_matrix)[2])){
    ests[[i]] = matrix(full_matrix[,i],dimension,dimension)
  }
  return(ests)
}
write_summary_measures = function(filename_ests,filename_error_measures,
                                  has_missingvalues=FALSE,Sigma){
  ests = read_ests(filename_ests)
  summary_measures = as.data.frame(sapply(1:length(ests),
                                          function(i) c(mean(abs(ests[[i]]-Sigma)),
                                                              sqrt(mean((ests[[i]]-Sigma)^2)))))
  rownames(summary_measures) = c("mean absolute error","rooted mean squared error")
  write.csv(summary_measures,file=filename_error_measures)
}

read_names_FITcomps_std_total = function(FITcomps_std_total, covar, model="no_missing_values"){
  all_min = rep(0,196)
  if(model=="no_missing_values"){
    all_min[which(!is.na(FITcomps_std_total[2,]))] = 1
    for(i in (2:12)){
      null_vec = rep(0,196)
      null_vec[which(!is.na(FITcomps_std_total[i,]))] = 1
      all_min = all_min*null_vec
    }
  } else{
    all_min[which(!is.na(FITcomps_std_total[12,]))] = 1
  }
  
  names_by_id = numeric(196)
  for(i in (1:(196*196))){
    if(names_by_id[covar$id_col[i]]!=0 && names_by_id[covar$id_col[i]]!=covar$name_o[i]){
      print(i)
      print("ERROR")
    }
    
    names_by_id[covar$id_col[i]] = covar$name_o[i]
  }

  return(list(names_by_id=names_by_id, all_min=all_min))
}

read_plot_FITcomps_std = function(filename){
  
  #Finding the matrices and epsilon
  dimension = 196
  num_observations = 12 # will be 11 later, because the first time interval is all missing values
  counter = matrix(0,nrow=num_observations,ncol=dimension)
  FITcomps_std_total = matrix(0,nrow=num_observations,ncol=dimension)
  #browser()
  
  n_cores=detectCores()
  
  
  res=mclapply(1:1000,function(i) {
    blabla = read.table(file=sprintf(filename,i), header=TRUE)
    blabla[is.na(blabla)]=0
    blabla
    },
               mc.cores = min(n_cores,8)) 
  
  FITcomps_std_total=Reduce("+",res)/length(res)
  FITcomps_std_total[FITcomps_std_total==0]=NA
  
  par(mfrow=c(2,6))
  all = matrix(NA,nrow=11,ncol=196)
  min = matrix(0,nrow=11,ncol=196)
  
  FITcomps_std = read.table(file=sprintf(filename,1), header= TRUE)
  
  for(i in (2:12)){
    print(paste0(rownames(FITcomps_std)[i]," ",sum(is.na(FITcomps_std_total[i,]))))
    hist(which(!is.na(FITcomps_std_total[i,])),breaks=(0:196)+exp(-16),main="",xlab=rownames(FITcomps_std)[i],ylab="not missing")
    all = c(all,which(!is.na(FITcomps_std_total[i,])))
  }
  
  par(mfrow=c(1,1))
  hist(all,breaks=(0:196)+exp(-16),main="total non-missing values by country",xlab="countries")
  print("percentage missing values:")
  print(mean(is.na(FITcomps_std_total[2:dim(FITcomps_std_total)[1],])))
  return(list(FITcomps_std_total=FITcomps_std_total,FITcomps_std=FITcomps_std))
}

preproc_FITcomps_std = function(all_min, names_by_id, FITcomps_std, covar){
  FITcomps_std_iso = numeric(length(which(all_min==1)))
  for(i in seq_along(FITcomps_std_iso)){
    FITcomps_std_iso[i] = as.numeric(substring(names(FITcomps_std)[which(all_min==1)][i],2,
                                               nchar(names(FITcomps_std)[which(all_min==1)][i])))
  }
  
  iso_id_key = numeric(196)
  for(i in (1:(196*196))){
    iso_id_key[covar$id_row[i]] = covar$iso_d[i]
  }
  
  adj_positions = numeric(length(FITcomps_std_iso))
  for(i in seq_along(FITcomps_std_iso)){
    adj_positions[i] = which(iso_id_key==FITcomps_std_iso[i])
  }
  
  names_by_id[adj_positions]
  
  
  dimension = dim(Gb)[1]
  num_observations = 12
  
  for(Fk in list(Gb,Gc)){
    cut_node = c()
    degree_Fk = rowSums(Fk*(diag(dimension)==0))
    for(i in (1:dimension)){
      counter=c()
      for(j in which((Fk*(diag(dimension)==0))[i,] == 1)){
        if(degree_Fk[i] != degree_Fk[j]){
          counter=c(counter,j)
        }
      }
      if(length(counter)>1){
        cut_node = c(cut_node,i)
      } else if(length(counter)==1){
        cut_node = counter[1]
      }
    }
    #browser()
    for(c_n in cut_node){
      adj_positions = adj_positions[which(adj_positions!=c_n)]
    } # remove cut-nodes
  }
  Fk = list()
  Fk[[1]] = Gb[adj_positions,adj_positions]
  diag(Gb)=1
  Fk[[2]] = Gc[adj_positions,adj_positions]
  Fk[[3]] = Gd[adj_positions,adj_positions]
  
  Gl = list()
  Ml = list()
  Al = list()
  vois <- sparseMatrix(i = covar$id_row, j = covar$id_col, x = covar$contig)
  
  
  M = diag(rowSums(vois))
  A = (diag(1/rowSums(vois))%*%vois)
  Gl[[1]] = matrix(0,ncol=196,nrow=196)
  Ml[[1]] = M
  Al[[1]] = A
  
  matList_final = list(Fk=Fk,Gl=Gl,Ml=Ml,Al=Al,vois=vois)
  return(list(matList_final=matList_final,adj_positions=adj_positions,iso_id_key=iso_id_key,FITcomps_std_iso=FITcomps_std_iso))
}