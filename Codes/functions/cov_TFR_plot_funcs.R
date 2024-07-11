PACKAGES_VIS = c("UpSetR", "ggplot2", "grid", "plyr", "ggheatmap",
                 "igraph", "ggcorrplot", "gridExtra", "corrplot", "RColorBrewer") # plot stuff
# PACKAGES_PARALLEL = c("foreach", 
#                       "parallel", "doParallel") # parallelization stuff
PACKAGES = c(PACKAGES_VIS) # all of the stuff
lapply(PACKAGES, require, character.only = TRUE)

plot_cov = function(matList,Sigma,SigmaHat_list,colvec,model="corY",
                    ests_names=NULL, order=1:length(ests_names)){
  
  covY = SigmaHat_list[[1]]
  #browser()
  p = dim(Sigma)[1]
  diag(Sigma)=NA
  #browser()
  colvec=colvec[order]
  ests_names=ests_names[order]
  
  data = cbind(Sigma[diag(p)!=1],sapply(seq_along(SigmaHat_list),function(i) c(SigmaHat_list[[i]][diag(p)!=1])))
  df = as.data.frame(data)
  #browser()
  
  # if(length(names(df))==5){
  #   names(df) = c("Sigma","IVE","Imputed","SCE")
  #   #browser()
  #   ggplot(df,aes(x=Sigma,y=Imputed,col="Imputed")) + geom_point() +
  #     geom_point(aes(x=Sigma,y=IVE,col="IVE"))+
  #     geom_point(aes(x=Sigma,y=SCE,col="SCE"))+
  #     scale_color_manual(values=c("black","brown","orange2"))+
  #     geom_line(aes(x=Sigma,y=Sigma),colour="darkgrey",linetype="solid", linewidth=.8)
  #     #+ theme(text=element_text(size=3))
  # } else{
  names(df) = c("Sigma",ests_names)#c("Sigma","Pearson","LW","Glasso","FM","hatSigma0","hatSigma")
  test2 = function(){
    for(i in seq_along(names(df)[-1])){
      df = cbind(df,i)
    }
    return(df)
  }
  df=test2()
  names(df) = c("Sigma",ests_names, paste0("V",1:(round((length(names(df))-1)/2))))
  df$Sigma2 = df$Sigma
    #browser()
  test = function(order){
    plot1 = ggplot(df,aes(x=Sigma,y=.data[[ests_names[[1]]]]))# + geom_point()
    for(i in order){
      plot1 = plot1 + geom_point(aes(x=Sigma,y=.data[[ests_names[[i]]]],col=as.factor(.data[[names(df)[length(ests_names)+1+i]]][1])))
    }
    return(plot1)
  }
  #colvec=c("brown","beige","pink","pink3","black","orange2","green")
  # scale_fill_manual(values = c("brown","grey","pink","pink3","beige","orange2","darkorange2"),
  #                   labels =  expression("Pearson","FM","Glasso","LW","IVE","SCE","WSCE"), 
  #                   name="")
  #browser()
  plot1 = test(order)+scale_color_manual(name="estimator",values=colvec,labels=ests_names)
  plot1 + geom_line(aes(x=Sigma2,y=Sigma2),colour="black",linetype="solid", linewidth=.8)+ylab("correlation")
}
  
  # plot(Sigma, covY, pch=20,xlab="Sigma",ylab=model); abline(a=0, b=1, h=0, v=0)
  # for(i in seq_along(SigmaHat_list)[-1]){
  #   points(Sigma, SigmaHat_list[[i]], pch=20, col=colvec[i-1])
  # }
  

#}

plot_cov_simple = function(id_min, matList, parm,corY,SigmaHat){
  #browser()
  matList_final = matList
  transformed_parm_01=sapply(1:length(c(parm)),function(i) c(parm)[[i]]) 
  SigmaHat1=SigmaHat
  
  diag(corY)=NA
  plot(corY,SigmaHat1)
  abline(h=transformed_parm_01[3],col="gold")
  abline(h=transformed_parm_01[2]+transformed_parm_01[3],col="green")
  
  rho = transformed_parm_01[5]
  mean_neighbor_effect = transformed_parm_01[6]#calc_mean_neighbor_effect(matList, rho, id_min)
  
  abline(h=mean_neighbor_effect+transformed_parm_01[2]+transformed_parm_01[3],col="purple")
  abline(h=mean_neighbor_effect+transformed_parm_01[2]+transformed_parm_01[3]+transformed_parm_01[1],col="red")
  
  
}

plot_param = function(matList, filename_true_param=NULL, filename_param_fit,
                      id_min = 1:dim(matList$Al[[1]])[1],error=TRUE,version="upsetplot"){
  #browser()
  parm_est = as.matrix(read_param(filename=filename_param_fit))
  
  n = dim(matList$Fk[[1]])[1]
  
  #use combined_effects if necessary
  if(length(parm_est)==8){
    Sigma_est <- CovMat_03(parm_est, matList,id_min=id_min, link=combined_matList)$Sigma
  } else{
    Sigma_est <- CovMat_03(parm_est, matList,id_min=id_min)$Sigma
  }
  
  #Figure 1: The different effects of sigma
  
  par(mfrow = c(1, 1))
  
  # Sigma_test = Sigma
  # diag(Sigma_test) = NA
  # plot(Sigma_test,cor(sim2$Y),xlab="Sigma",ylab="corY")
  # points(Sigma_test[matList2$Fk[[1]]==1],cor(sim2$Y)[matList2$Fk[[1]]==1],col="RED")
  # points(Sigma_test[matList2$Fk[[2]]==1],cor(sim2$Y)[matList2$Fk[[2]]==1],col="PURPLE")
  # points(Sigma_test[matList2$Al[[1]]!=0],cor(sim2$Y)[matList2$Al[[1]]!=0],col="GREEN")
  
  if(!is.null(filename_true_param)){
    parm = as.matrix(read_param(filename=filename_true_param))
    #browser()
    if(error){
      Sigma <- CovMat_03(parm, matList,id_min=id_min, link=combined_matList)$Sigma
      df = as.data.frame(matrix(ncol=6,nrow=n*(n-1)/2))
      names(df) = c("error","correlation","comcol","region","neighbor","nothing")
      #browser()
      df$error = c(Sigma[upper.tri(Sigma)]-Sigma_est[upper.tri(Sigma_est)])
    } else{
      parm=c(.05,.09,.11,.26)
      Sigma <- CovMat_03(parm, matList,id_min=id_min, combined_effects="FosdickRaftery")$Sigma
      df = as.data.frame(matrix(ncol=6,nrow=n*(n-1)/2))
      #browser()
      names(df) = c("estimate","correlation","comcol","region","neighbor","nothing")
      #browser()
      df$estimate = c(Sigma_est[upper.tri(Sigma_est)])
    }
    #browser()
    df$correlation = c(Sigma[upper.tri(Sigma)])
    df$comcol = c(matList$Fk[[1]][upper.tri(Sigma)])
    df$region = c(matList$Fk[[2]][upper.tri(Sigma)])
    matList$Al[[1]][is.na(matList$Al[[1]])]=0
    df$neighbor = c(((matList$Al[[1]][id_min,id_min]!=0)+0)[upper.tri(Sigma)])
    df$nothing = (1-df$comcol)*(1-df$region)*(1-df$neighbor)
    if(error){
      #browser()
      upset(df, boxplot.summary = c("correlation","error"))
    } else{
      #browser()
      #plot(Sigma_est,Sigma)
      if(version=="scatterplot"){
        plot(df$estimate,df$correlation)
        points(df$estimate[df$comcol==1],df$correlation[df$comcol==1],col="red")
        points(df$estimate[df$region==1],df$correlation[df$region==1],col="green3")
        points(df$estimate[df$neighbor==1],df$correlation[df$neighbor==1],col="purple")
        abline(a=0,b=1,col="black")
      } else{
        upset(df, boxplot.summary = c("correlation","estimate"))
      }
      #points(df$estimate[df$region],df$correlation[df$region],col="green")
    }
  } else{
    df = as.data.frame(matrix(ncol=5,nrow=n*(n-1)/2))
    names(df) = c("correlation","comcol","region","neighbor","nothing")
    df$correlation = c(Sigma_est[upper.tri(Sigma_est)])
    df$comcol = c(matList$Fk[[1]][upper.tri(Sigma_est)])
    df$region = c(matList$Fk[[2]][upper.tri(Sigma_est)])
    matList$Al[[1]][is.na(matList$Al[[1]])]=0
    df$neighbor = c(((matList$Al[[1]][id_min,id_min]!=0)+0)[upper.tri(Sigma_est)])
    df$nothing = (1-df$comcol)*(1-df$region)*(1-df$neighbor)
    #browser()
    upset(df, boxplot.summary = c("correlation"))
  }
  #browser()
  
  #library(devtools)
  #install_github("jokergoo/ComplexHeatmap")
  #library(ComplexHeatmap)
  # So far: No column where variables have "no effect"; We might want that?
}

plot_data = function(matList, names_by_id, id_min, show_names=TRUE){
  
  #browser()
  par(mfrow=c(1,1))
  for(k in (1:2)){
    diag(matList$Fk[[k]])=0
    g3 <- graph_from_adjacency_matrix(as.matrix(matList$Fk[[k]]),mode="undirected")
    #browser()
    V(g3)$name <- names_by_id[id_min]#substring(names_by_id[id_min],1,2)
    if(show_names){
      plot(g3, vertex.size =.1,vertex.label.cex=3/5)
    } else{
      plot(g3, vertex.size =.1,vertex.label.cex=3/15)
    }
    diag(matList$Fk[[k]])=1
  }
  
  par(mfrow=c(1,1))
  A = matList$Al[[1]][id_min,id_min]
  A[is.na(A)] = 0
  g3 <- graph_from_adjacency_matrix(as.matrix((A!=0)+0),mode="undirected")
  V(g3)$name <- names_by_id[id_min]
  if(show_names){
    plot(g3,vertex.size=.1,vertex.label.cex=3/5)
  } else{
    #browser()
    plot(g3,vertex.size=.1,vertex.label.cex=3/15)
  }
  
  # only neighbor, but not comcol or region
  par(mfrow=c(1,1))
  g3 <- graph_from_adjacency_matrix(as.matrix((A!=0)*(matList$Fk[[1]]==0)*(matList$Fk[[2]]==0)+0),mode="undirected")
  V(g3)$name <- names_by_id[id_min]
  if(show_names){
    plot(g3,vertex.size=.1,vertex.label.cex=3/5)
  } else{
    #browser()
    plot(g3,vertex.size=.1,vertex.label.cex=3/15)
  }
  
  # neighbor + comcol
  par(mfrow=c(1,1))
  amatrix = as.matrix((A!=0)*(matList$Fk[[1]]==1))
  no_islands_id = apply(amatrix,1, function(a) sum(a)!=0)
  islands_id = (1-no_islands_id)==1
  g3 <- graph_from_adjacency_matrix(amatrix,mode="undirected")
  V(g3)$name <- names_by_id[id_min]
  V(g3)$name[islands_id] = ""
  plot(g3,vertex.size=.1,vertex.label.cex=3/15)

  # neighbor + region
  par(mfrow=c(1,1))
  amatrix = as.matrix((A!=0)*(matList$Fk[[2]]==1))
  no_islands_id = apply(amatrix,1, function(a) sum(a)!=0)
  islands_id = (1-no_islands_id)==1
  g3 <- graph_from_adjacency_matrix(amatrix,mode="undirected")
  V(g3)$name <- names_by_id[id_min]
  V(g3)$name[islands_id] = ""
  plot(g3,vertex.size=.1,vertex.label.cex=3/15)
  
  # neighbor + comcol - region
  par(mfrow=c(1,1))
  amatrix = as.matrix((A!=0)*(matList$Fk[[1]]==1)*(matList$Fk[[2]]==0))
  no_islands_id = apply(amatrix,1, function(a) sum(a)!=0)
  islands_id = (1-no_islands_id)==1
  g3 <- graph_from_adjacency_matrix(amatrix,mode="undirected")
  V(g3)$name <- names_by_id[id_min]
  V(g3)$name[islands_id] = ""
  plot(g3,vertex.size=.1,vertex.label.cex=3/15)
  
  
  
}

plot_sims = function(sims_errors_and_bic,filename, has_missingvalues=FALSE){
  browser()
  par(mfrow=c(2,1))
  # if(has_missingvalues){
  #   mat <- cbind(sims_errors_and_bic$mae1.2,
  #                sims_errors_and_bic$mae1.1,
  #                sims_errors_and_bic$mae1.3)
  #   names = c("Imputed","IVE","SCE") 
  #   mat2 <- cbind(sims_errors_and_bic$rmse1.2,
  #                 sims_errors_and_bic$rmse1.1,
  #                 sims_errors_and_bic$rmse1.3)
  mat <- cbind(sims_errors_and_bic$mae1.1,
               sims_errors_and_bic$mae1.2,
               sims_errors_and_bic$mae1.3,
               sims_errors_and_bic$mae1.4,
               sims_errors_and_bic$mae1.5,
               sims_errors_and_bic$mae1.6,
               sims_errors_and_bic$mae1.7)
  names = ESTS_NAMES    

  #browser()
  colnames(mat) = names
  df1=as.data.frame(mat)
  library(reshape2)
  df1=melt(df1)
  df1$MAE = df1$value
  df1$estimator=df1$variable
  plot1 =ggplot(df1,aes(x=estimator,y=MAE)) + geom_boxplot()
  #boxplot(mat,main="Mean Absolute Error (MAE)")
  # colnames(mat2) = names
  # df2 = as.data.frame(mat2)
  # df2=melt(df2)
  # df2$RMSE=df2$value
  # df2$estimator = df2$variable
  # plot2 = ggplot(df2,aes(x=estimator,y=RMSE)) + geom_boxplot()
  
  #plot =  arrangeGrob(plot1,plot2,ncol=1)
  #ggsave(plot,filename=filename, width=5.3,height=4.07,device="pdf")
  #browser()
  ggsave(plot1,filename=filename, width=5.3,height=4.07,device="pdf")
  #boxplot(mat2,main="Rooted Mean Squared Error (RMSE)")
}

plot_heatmaps = function(matList, 
                         Sigma, Sigma2=NULL, Sigma3=NULL, 
                         filename, show_regions=FALSE,
                         names_by_id=NULL, id_min=NULL,
                         iso3=NULL, info_coutries=NULL, Sigma_true=Sigma,
                         Sigma_cluster=Sigma, main_range=NULL, main_color=NULL){
  #browser()
  set.seed(1)
  dummy_Sigma = data.frame(Sigma_true)
  names(dummy_Sigma)=(1:dim(Sigma_true)[1])
  
  hclust_cluster = hclust(as.dist(1-Sigma_true))
  
  ### TEST ###
  
  library(mclust)
  
  plotCAH <- function(res) {
    
    W = c(0, res$height)
    taille = length(W)
    plot(1:taille, W[taille:1], type="b", col="red")
  }
  
  # first method ! 
  resCAH = hclust_cluster # CAH avec Ward
  
  plot(resCAH)
  
  plotCAH(resCAH)
  K1 = 5 # number of clusters (from visualisation)
  abline(v=K1)

  cl1 = cutree(resCAH, K1)
  #browser()
  
  #adjustedRandIndex(cl1, cl2)
  
  ### TEST ###
  h_order = hclust(as.dist(1-Sigma_true))$order
  order_x = h_order
  order_y = h_order
  
  diag(Sigma)=0
  #browser()
  if(show_regions==TRUE){
    
    browser()
    #matList$Al[[1]][is.na(matList$Al[[1]])]=0
    #Sigma = ceiling(as.matrix(Sigma*matList$Al[[1]][id_min,id_min]))
    df = Sigma[unique(order_y),unique(order_y)][1:length(unique(order_y)),length(unique(order_y)):1]
    browser()
    df = df 
    
    df = as.data.frame(df)
    #df = Sigma_countries
    
    #rownames(df) = paste("s ",1:dim(Sigma)[1],sep = "")
    #colnames(df) = paste("s ",1:dim(Sigma)[1],sep = "")
    rownames(df)=paste(dim(Sigma)[1]:1,(names_by_id[id_min])[unique(order_y)])#[length(unique(order_y)):1])
    colnames(df)=paste(1:dim(Sigma)[1],(names_by_id[id_min])[unique(order_y)])#[length(unique(order_y)):1])
    
    browser()
    #comcol_mat = cbind(1-(diag(matList$Ml[[1]])>min(diag(matList$Ml[[1]]))),diag(matList$Ml[[1]])>min(diag(matList$Ml[[1]])))
    #comcol_mat = comcol_mat %*% t(comcol_mat)
    #comcol_mat = comcol_mat[id_min,id_min][unique(order_y),unique(order_y)]
    comcol_mat = matList$Fk[[1]][unique(order_y),unique(order_y)]
    # 
    # # find the different colonizers and regions
    counter = 1
    comcol_vec = c(1)
    comcol_res_vec = numeric(dim(Sigma)[1])
    comcol_res_vec[comcol_mat[1,]==1] = as.character(counter)#numeric(dim(Sigma)[1])
     
    for(i in 1:dim(Sigma)[1]){
      if(sum(sapply(comcol_vec,function(s) s%in%which(comcol_mat[i,]==1)))==0){
        counter = counter + 1
        comcol_vec[counter] = i
        comcol_res_vec[comcol_mat[i,]==1]=as.character(counter)
      }
    }

    # reorder
    comcol_res_vec = sapply(seq_along(comcol_res_vec),
                            function(s) as.character(which(unique(comcol_res_vec)==comcol_res_vec[s])[1]))
    
    
    ### TEST ###
    
    # reg_mat = matList$Fk[[2]][unique(order_y),unique(order_y)]
    # # 
    # # # find the different colonizers and regions
    # counter = 1
    # reg_vec = c(1)
    # reg_res_vec = numeric(dim(Sigma)[1])
    # reg_res_vec[reg_mat[1,]==1] = as.character(counter)#numeric(dim(Sigma)[1])
    # 
    # for(i in 1:dim(Sigma)[1]){
    #   if(sum(sapply(reg_vec,function(s) s%in%which(reg_mat[i,]==1)))==0){
    #     counter = counter + 1
    #     reg_vec[counter] = i
    #     reg_res_vec[reg_mat[i,]==1]=as.character(counter)
    #   }
    # }
    # 
    # # reorder
    # reg_res_vec = sapply(seq_along(reg_res_vec),
    #                         function(s) as.character(which(unique(reg_res_vec)==reg_res_vec[s])[1]))
    #browser()
    
    
    ### TEST ###
    
    browser()
    
    id_nocol = as.numeric(comcol_res_vec[sapply(names(df), function(s) grepl("France",s))])
    id_fra = as.numeric(comcol_res_vec[sapply(names(df), function(s) grepl("Senegal",s))])
    id_usa = as.numeric(comcol_res_vec[sapply(names(df), function(s) grepl("Philippines",s))])
    id_gbr = as.numeric(comcol_res_vec[sapply(names(df), function(s) grepl("India",s))])
    id_prt = as.numeric(comcol_res_vec[sapply(names(df), function(s) grepl("Cape Verde",s))])
    id_rus = as.numeric(comcol_res_vec[sapply(names(df), function(s) grepl("Estonia",s))])
    id_nld = as.numeric(comcol_res_vec[sapply(names(df), function(s) grepl("Indonesia",s))])
    id_esp = as.numeric(comcol_res_vec[sapply(names(df), function(s) grepl("Western Sahara",s))])
    id_bel = as.numeric(comcol_res_vec[sapply(names(df), function(s) grepl("Burundi",s))])
    id_jpn = as.numeric(comcol_res_vec[sapply(names(df), function(s) grepl("Republic of Korea",s))])
    id_jpn = id_jpn[id_jpn!=id_nocol] # Republic of Korea is in two names


    letters[id_nocol] = "."
    letters[id_gbr] ="GBR"
    letters[id_nld] = "NLD"
    letters[id_usa] = "USA"
    letters[id_prt] = "PRT"
    letters[id_fra] = "FRA"
    letters[id_esp] = "ESP"
    letters[id_bel] = "BEL"
    letters[id_rus] = "RUS"
    letters[id_jpn] = "JPN"

    comcol_res_vec=sapply(as.numeric(comcol_res_vec),function(s) letters[s])
    #browser()

    Sigma_countries <- Sigma_true %>% 
      as.data.frame() %>% 
      mutate(iso3 =  iso3[id_min] ) %>% 
      left_join(info_coutries, by =c( "iso3" = "ISO3 Alpha-code"))
    #sapply(1:K1,function(s) length(Sigma_countries$country_name[cl1==s]))
    
    continent_names = rep("others",195)
    #list_1 = c("Northern Europe","Western Africa","Eastern Asia", "Western Asia", "Central America", "Polynesia", "Australia/New Zealand", "Northern Africa", "Northern America", "Middle Africa")
    #africa_coords = which(sapply(1:195, function(s) Sigma_countries$reg_name[s]%in%list_1))
    africa_coords = which(sapply(1:195, function(s) grepl("Africa",Sigma_countries$reg_name[s])))
    continent_names[africa_coords] = "Africa"
    #continent_names[africa_coords] = "Block1"
    #list_2 = c("Caribbean","Eastern Africa","Western Europe","Polynesia", "Southern Europe", "Melanesia", "Micronesia", "Northern Asia","Southern Africa","Western Asia")
    asia_coords = which(sapply(1:195, function(s) grepl("Asia",Sigma_countries$reg_name[s])))
    continent_names[asia_coords] = "Asia"
    #asia_coords = which(sapply(1:195, function(s) Sigma_countries$reg_name[s]%in%list_2))
    #continent_names[asia_coords] = "Block2"
    
    #list_3 = c("Eastern Europe", "South-Eastern Asia", "South America")
    america_coords = which(sapply(1:195, function(s) grepl("America",Sigma_countries$reg_name[s])))
    continent_names[america_coords] = "America"
    #america_coords = which(sapply(1:195, function(s) Sigma_countries$reg_name[s]%in%list_3))
    #continent_names[america_coords] = "Block3"
    #list_4 = c("South-Central Asia")
    europe_coords = which(sapply(1:195, function(s) grepl("Europe",Sigma_countries$reg_name[s])))
    continent_names[europe_coords] = "Europe"
    #europe_coords = which(sapply(1:195, function(s) Sigma_countries$reg_name[s]%in%list_4))
    #continent_names[europe_coords] = "Block4"
    continent_names_old = continent_names
    cluster_names = c("",cl1)[-1]
    
    browser()
    
    cluster_res_vec = cluster_names[unique(order_y)]
    continent_res_vec = continent_names[unique(order_y)]
    reg_res_vec=Sigma_countries$reg_name[unique(order_y)]#[length(unique(order_y)):1]
    comcol_res_vec=comcol_res_vec[length(comcol_res_vec):1]#[unique(order_y)]#[length(unique(order_y)):1]
    #comcol_res_vec=Sigma_countries$short_colonizer1[unique(order_y)][length(unique(order_y)):1]
    
    row_metaData = data.frame(comcols=comcol_res_vec, clusters=cluster_res_vec[195:1])
    col_metaData = data.frame(regions=reg_res_vec, continents=continent_res_vec)
    
    regcol <- turbo(length(unique(reg_res_vec)))
    names(regcol) <- unique(reg_res_vec)
    
    comcolcol <- turbo(length(unique(comcol_res_vec)))#c("#98D352","#FF7F0E")
    names(comcolcol) <- unique(comcol_res_vec)#c("Normal","Tumor")
    
    clustercol = turbo(length(unique(cluster_res_vec)))
    names(clustercol) = unique(cluster_res_vec)
    
    continentcol = turbo(length(unique(continent_res_vec)))
    names(continentcol) = unique(continent_res_vec)
    
    col <- list(comcols=comcolcol,
                regions=regcol,
                clusters=clustercol,
                continents=continentcol)
    
    # add zeroes to keep the order
    df_col = colnames(df)#[h_order]
    browser()
    df_col[as.numeric(sapply(df_col,function(s) str_split(s,pattern=" ")[[1]][1]))<10] = 
      paste("00",df_col[as.numeric(sapply(df_col,function(s) str_split(s,pattern=" ")[[1]][1]))<10],sep="")
    df_col[(as.numeric(sapply(df_col,function(s) str_split(s,pattern=" ")[[1]][1]))>=10) & 
             (as.numeric(sapply(df_col,function(s) str_split(s,pattern=" ")[[1]][1]))<100)] =
      paste("0",    df_col[(as.numeric(sapply(df_col,function(s) str_split(s,pattern=" ")[[1]][1]))>=10) & 
                             (as.numeric(sapply(df_col,function(s) str_split(s,pattern=" ")[[1]][1]))<100)],sep="")

    
    #rownames(df)=colnames(df)
    rownames(df) <- colnames(df) <- df_col[length(df_col):1] # flip y-axis
    rownames(row_metaData) <- rownames(df)
    rownames(col_metaData) <- colnames(df)
    
    
    
    browser()
    
    ### TEST ###
    h_order = hclust(as.dist(1-Sigma_cluster))$order
    df=as.data.frame(Sigma)
    # browser()
    names(df) = names_by_id[id_min]
    df = df[h_order, h_order][195:1,]
    # rownames(df) <- colnames(df) <- rownames(row_metaData)#[length(rownames(row_metaData)):1]
    # 
    # 
    # #row_metaData = row_metaData[h_order,]
    # #col_metaData = col_metaData[h_order,]
    # 
    df_col = paste0(195:1," ",colnames(df))#[h_order]
    # browser()
    df_col[as.numeric(sapply(df_col,function(s) str_split(s,pattern=" ")[[1]][1]))<10] =
      paste("00",df_col[as.numeric(sapply(df_col,function(s) str_split(s,pattern=" ")[[1]][1]))<10],sep="")
    df_col[(as.numeric(sapply(df_col,function(s) str_split(s,pattern=" ")[[1]][1]))>=10) &
             (as.numeric(sapply(df_col,function(s) str_split(s,pattern=" ")[[1]][1]))<100)] =
      paste("0",    df_col[(as.numeric(sapply(df_col,function(s) str_split(s,pattern=" ")[[1]][1]))>=10) &
                             (as.numeric(sapply(df_col,function(s) str_split(s,pattern=" ")[[1]][1]))<100)],sep="")
    # 
    # 
    # #rownames(df)=colnames(df)
    rownames(df) <- colnames(df) <- df_col[length(df_col):1] # flip y-axis
    # rownames(row_metaData) <- rownames(df)
    # rownames(col_metaData) <- colnames(df)
    
    browser()
    row_metaData = row_metaData[sort(sapply(rownames(row_metaData), function(s) str_split(s," ")[[1]][2]),index.return=TRUE)$ix,][order(sort(sapply(rownames(df), function(s) str_split(s," ")[[1]][2]),index.return=TRUE)$ix),]
    #row_metaData = row_metaData[195:1,]
    
    col_metaData = as.data.frame(cbind(Sigma_countries$reg_name,Sigma_countries$reg_name))
    rownames(col_metaData) = Sigma_countries$country_name
    col_metaData = col_metaData[sort(rownames(col_metaData),index.return=T)$ix,][order(sort(sapply(colnames(df), function(s) str_split(s," ")[[1]][2]),index.return=TRUE)$ix),]
    colnames(col_metaData)=c("regions","continents")
    
    rownames(row_metaData) <- rownames(df)
    rownames(col_metaData) <- colnames(df)
    
    Sigma_countries$country_name[Sigma_countries$reg_name=="Western Europe"]
    
    ### TEST ###
    browser()

    col_metaData$regions = col_metaData$regions[195:1]
    row_metaData$clusters=NULL
    col_metaData$continents=NULL
    #rownames(df)<- colnames(df)<-rownames(row_metaData)<-rownames(col_metaData)<-rownamesdf[195:1]#[h_order]
    if(!is.null(main_range)){
      df[1,195]=main_range[1]
      df[2,194]=main_range[2]
    }
    #df[1,195]=1
    #df[1,1]=-1
    if(sum(df<0)!=0){
      browser()
      if(is.null(main_color)){
        main_color=c(mako(2*sum(df<0), alpha = 1, begin = 0, end = 1, direction = -1)[(1+sum(df<0)):(2*sum(df<0))],
                     magma(sum(df>=0), alpha = 1, begin = 0, end = 1, direction = 1))
      }

      p = ggheatmap::ggheatmap(df, annotation_rows = row_metaData,
                               annotation_cols = col_metaData,
                               annotation_color = col,
                               #text_show_rows = c(""),
                               text_show_cols = c(""),
                               legendName="corr",
                               color=main_color)
    } else{
      if(is.null(main_color)){
        main_color = magma(sum(c(Sigma)>=0), alpha = 1, begin = 0, end = 1, direction = 1)
      }

      p = ggheatmap::ggheatmap(df, annotation_rows = row_metaData,
                               annotation_cols = col_metaData,
                               annotation_color = col,
                               #text_show_rows = c(""),
                               text_show_cols = c(""),
                               legendName="corr",
                               color=main_color)
    }

    
    # p = ggheatmap::ggheatmap(df, annotation_rows = row_metaData,
    #               annotation_cols = col_metaData,
    #               annotation_color = col,
    #               #text_show_rows = c(""),
    #               text_show_cols = c(""),
    #               legendName="corr",
    #               color=c(cividis(sum(c(Sigma)>=0), alpha = 1, begin = 0, end = 1, direction = -1)[sum(c(Sigma)>=0):(sum(c(Sigma)>=0)-sum(c(Sigma)<0)+398)],
    #                       magma(sum(c(Sigma)>=0), alpha = 1, begin = 0, end = 1, direction = 1)[1:(sum(c(Sigma)>=0))]))
    #browser()
    ggsave(filename, p, device="jpeg", width=20, height=20)
    
    library(cowplot)
    legend1 <- cowplot::get_legend(p[[1]])
    legend2 <- cowplot::get_legend(p[[2]])
    legend3 <- cowplot::get_legend(p[[3]])
    plot2=grid.arrange(legend2[[1]][[1]],legend3[[1]][[1]],legend1[[1]][[1]])
    ggsave("atelier/sim_final_n195_full_data_covmat_legend.jpeg",plot2, device="jpeg", width=4, height=9.5)

    browser()
    print(median(df[df<0]))
    
    return(list(main_range=range(df),main_color=main_color))
    
  } else{#plot 3 matrices (Adrian's, ours without-, ours with interaction effects)
    if(!is.null(Sigma2)){
      diag(Sigma2)=0
      diag(Sigma3)=0
      plotSigma = ggcorrplot(Sigma[unique(order_y),unique(order_y)], outline.col = "white") + 
        scale_fill_viridis(option="magma",limit=c(0,max(Sigma,Sigma2,Sigma3))) + ggtitle("Model from Fosdick,Raftery")
      plotSigma2 = ggcorrplot(Sigma2[unique(order_y),unique(order_y)], outline.col = "white") +
        scale_fill_viridis(option="magma",limits=c(min(Sigma),max(Sigma,Sigma2,Sigma3))) + ggtitle("Base Model")
      plotSigma3 = ggcorrplot(Sigma3[unique(order_y),unique(order_y)], outline.col = "white") +
        scale_fill_viridis(option="magma",limits=c(min(Sigma),max(Sigma,Sigma2,Sigma3))) + ggtitle("Interaction Effects Model")
      
      #browser()
      
      if(dim(Sigma)[1]==195){#needs more space if it is large
        #browser()
        ggsave(filename,
               grid.arrange(plotSigma, plotSigma2, plotSigma3,ncol=2),
               device = "jpeg",width=15.6,height=12.14)
        return(cl1)
      } else{
        #browser()
        ggsave(filename,
               grid.arrange(plotSigma, plotSigma2, plotSigma3,ncol=3),
               device = "jpeg",width=10.6,height=8.14)
      }
      
    } else{
      plotSigma = ggcorrplot(Sigma[unique(order_y),unique(order_y)][1:length(unique(order_y)),length(unique(order_y)):1],
                             outline.col = "white",lab_size = 0) + 
        scale_fill_viridis(option="magma") + ggtitle("Correlation Matrix") + 
        theme(axis.text.y=element_blank(),axis.text.x=element_blank())
      plotF1 = ggcorrplot(as.matrix(matList$Fk[[1]])[unique(order_y),unique(order_y)][1:length(unique(order_y)),length(unique(order_y)):1]
                          , outline.col = "white") +
        scale_fill_viridis(option="magma") + ggtitle("Common Colonizers") +
        theme(axis.text.y=element_blank(),axis.text.x=element_blank())
      plotF2 = ggcorrplot(as.matrix(matList$Fk[[2]])[unique(order_y),unique(order_y)][1:length(unique(order_y)),length(unique(order_y)):1]
                          , outline.col = "white") +
        scale_fill_viridis(option="magma") + ggtitle("Regions") + 
        theme(axis.text.y=element_blank(),axis.text.x=element_blank())
      plotAl = ggcorrplot(((as.matrix(matList$Al[[1]])!=0)+0)[unique(order_y),unique(order_y)][1:length(unique(order_y)),length(unique(order_y)):1],
                          outline.col = "white") +
        scale_fill_viridis(option="magma") + ggtitle("Neighbors") +
        theme(axis.text.y=element_blank(),axis.text.x=element_blank())
      
      
      
      ggsave(filename,
             grid.arrange(plotSigma, plotF1, plotF2, plotAl,ncol=2),
             device = "jpeg",width=10.6,height=8.14)#dpi=700
    }
  }
  
}