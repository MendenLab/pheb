library(stringr)

# Convert a beta matrix to m-values
beta2m <- function(x){
    return(log2(x/(1- x)))
}

# Load R data object and immediately return it
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

sd_BT <- function(vector){
 output <- lapply(1:length(vector), function(x) sd(vector[-x],na.rm = T))
 return(min(unlist(output)))
}


cut_df <- function(
  ### cuts a dataframe in columns and returns the one with a certain index
  ######################################################  
  df = NULL,
  totaln = NULL,
  index = NULL
){
  tmp <- cut(1:ncol(df), breaks = totaln, labels = F)
  df <- df[,tmp == index, drop = F]
  return(df)
}######################################################


enrichment <- function(
  ###
  ### Form an abundance table, calculate enrichment with hypergeometric test
  ### 
  ######################################################
  tab, # table of counts for two covariates
  col_label, # name of columns
  row_label, # name of rows
  digits # round of p-value to digit x
){
  tab_test <- matrix(ncol = ncol(tab), nrow = nrow(tab))
  for(i in 1:nrow(tab_test)){
    for(j in 1:ncol(tab_test)){
      p_under <- phyper(q = tab[i,j],
                        m = sum(tab[i,]),
                        n = sum(tab[-i,]),
                        k = sum(tab[,j]), 
                        lower.tail = T # Under-representation
      )
      p_over <- phyper(q = tab[i,j]-1,
                       m = sum(tab[i,]),
                       n = sum(tab[-i,]),
                       k = sum(tab[,j]), 
                       lower.tail = F # Over-representation
      )
      p <- which.min(c(p_under,p_over))
      if(p==1){p <- -p_under}
      if(p==2){p <- p_over}
      tab_test[i,j] <-p 
    }
    row.names(tab_test) = row.names(tab)
    colnames(tab_test)=colnames(tab)
  }
  tab_test <- reshape2::melt(tab_test); colnames(tab_test) = c(row_label,col_label,"value")
  tab_test$label <- abs(round(tab_test$value, digits =digits))
  tab_test$label <- unlist(lapply(tab_test$label, function(x) if(x < 0.05){paste0("p=",as.character(format.pval(x)),"*")}else{paste0("p=",as.character(format.pval(x)),"")}))
  tab_test$Enrichment <- sign(tab_test$value)%>%as.character%>%dplyr::recode("1"="Over-enrichment","-1"="Under-enrichment")
  tab_test$p <- abs(tab_test$value)
  tab_test$fdr <- p.adjust(abs(tab_test$value), method = "BH")
  tab_test$abundance <- lapply(1:nrow(tab_test), function(x) 
    list(q = tab[tab_test[x,1],tab_test[x,2]]-1,m = sum(tab[tab_test[x,1],]),n = sum(tab[-c(tab_test[x,1]),]),k = sum(tab[,tab_test[x,2]]) )
  )
  
  
  return(tab_test)
}


#'@return the shortest path as a list of vertices or NULL if there is no path between src and dest
shortest_path <- function(graph, src, dest){
  path <- suppressWarnings(get.shortest.paths(graph, src, dest))
  path <- names(path$vpath[[1]])
  if (length(path)==1) NULL else path
} 

#'@return the sum of the weights of all the edges in the given path
path_weight <- function(path, graph) sum(E(graph, path=path)$weight)

#'@description sorts a list of paths based on the weight of the path
sort_paths <- function(graph, paths) paths[paths %>% sapply(path_weight, graph) %>% order]

#'@description creates a list of edges that should be deleted
find_edges_to_delete <- function(A,i,rootPath){
  edgesToDelete <- NULL
  for (p in A){
    rootPath_p <- p[1:i]
    if (all(rootPath_p == rootPath)){
      edge <- paste(p[i], ifelse(is.na(p[i+1]),p[i],p[i+1]), sep = '|')
      edgesToDelete[length(edgesToDelete)+1] <- edge
    }
  }
  unique(edgesToDelete)
}

#returns the k shortest path from src to dest
#sometimes it will return less than k shortest paths. This occurs when the max possible number of paths are less than k
k_shortest_yen <- function(graph, src, dest, k){
  if (src == dest) stop('src and dest can not be the same (currently)')
  
  #accepted paths
  A <- list(shortest_path(graph, src, dest))
  if (k == 1) return (A)
  #potential paths
  B <- list()
  
  for (k_i in 2:k){
    prev_path <- A[[k_i-1]]
    num_nodes_to_loop <- length(prev_path)-1
    for(i in 1:num_nodes_to_loop){
      spurNode <- prev_path[i]
      rootPath <- prev_path[1:i]
      
      edgesToDelete <- find_edges_to_delete(A, i,rootPath)
      t_g <- delete.edges(graph, edgesToDelete)
      #for (edge in edgesToDelete) t_g <- delete.edges(t_g, edge)
      
      spurPath <- shortest_path(t_g,spurNode, dest)
      
      if (!is.null(spurPath)){
        total_path <- list(c(rootPath[-i], spurPath))
        if (!total_path %in% B) B[length(B)+1] <- total_path
      }
    }
    if (length(B) == 0) break
    B <- sort_paths(graph, B)
    A[k_i] <- B[1]
    B <- B[-1]
  }
  A
}


### Function: Plot association
plot_asso <- function(df = df, # data.frame
                      meth_df = results_collected, 
                      gex_df = cancertypes, 
                      cpg_column = "cpg",
                      gene_column = "gdsc_eqtl",
                      drug_column = "Drug",
                      cancertype_column = "cancertype",
                      drug_column_label = "drug",
                      index_meth = 1,
                      index_gex = 1,
                      drug_info = drugs_info,
                      region_column = "elmer_region",
                      gex_y = TRUE, #if FALSE then response is on y
                      use_tcga = F
){ 
  
  #results_collected mine here for getting resp
  to <- list()
  p <- list()
  for(i in 1:nrow(df)){
    meth <- meth_df[[df[i,cancertype_column]]]
    if(!use_tcga){
      resp <- meth$DAT$resp
      meth <- meth$DAT$meth
    }
    
    if(is.null(meth)){stop("methylation dataframe is empty !")}
    if(!is.list(gex_df)){
      gex <- gex_df[intersect(row.names(meth),row.names(gex_df)),]
      meth <- meth[intersect(row.names(meth),row.names(gex)),]
    }else{
      gex <- gex_df[[df$cancertype[i]]]
      
      # adjust TCGA rownames
      if(use_tcga){
        tmp <- unlist(lapply(row.names(gex), function(x) paste(strsplit(x,"-")[[1]][1:4],collapse="-")))
        tmp.duplicated <- duplicated(tmp)
        gex <- gex[!tmp.duplicated,]
        row.names(gex) = tmp[!tmp.duplicated]
        
        tmp <- unlist(lapply(row.names(meth), function(x) paste(strsplit(x,"-")[[1]][1:4],collapse="-")))
        tmp.duplicated <- duplicated(tmp)
        meth <- meth[!tmp.duplicated,]
        row.names(meth) = tmp[!tmp.duplicated]
        
      }
      ######################
      gex <- gex[intersect(row.names(meth),row.names(gex)),]
      meth <- meth[intersect(row.names(meth),row.names(gex)),]
    }
    
    #sort 
    resp <- tryCatch(resp[intersect(row.names(meth),row.names(resp)),], error = function(e) NULL)
    if(!use_tcga){
      gex_new <- tryCatch(gex[intersect(row.names(resp),row.names(gex)),], error = function(e) NULL)
    }else{
      gex_new <- gex #full
    }
    if(!is.null(gex_new)){gex <- gex_new}
    if(!use_tcga){
      meth_new <- tryCatch(meth[intersect(row.names(resp),row.names(meth)),], error = function(e) NULL)
    }else{
      meth_new <- meth #full
    }
    if(!is.null(meth_new)){meth <- meth_new}
    
    #extract info
    ct <- df[i, cancertype_column]
    message(paste0("Cancer type: ",ct))
    drug <- df[i,drug_column]
    print(drug)
    drug_label <- df[i,drug_column_label]
    gene <- df[i,gene_column]
    message(paste0("Total number of genes: ",length(unlist(gene))))
    gene <- unlist(gene)[min(index_gex,length(unlist(gene)))]
    site <- unlist(df[i, cpg_column])
    site <- site[site %in% colnames(meth)]
    message(paste0("Total number of sites: ",length(unlist(site))))
    
    if(all(index_meth == "all")){
      message("Getting all cpgs..")
      site <- unlist(site)
    }else{
      site <- unlist(site)[min(index_meth, length(unlist(site)))]
    }
    
    resp <- tryCatch(resp[,drug, drop = F], error = function(e) NULL)
    gex <- gex[,intersect(gene, colnames(gex)), drop = F]
    meth <- meth[,site, drop = F]
    
    print(dim(gex))
    print(dim(meth))
    print(dim(resp))
    region <- df[i,region_column]
    region <- gsub("-","",region)
    if(use_tcga){
      to[[i]] <- cbind(meth, gex)
    }else{
      to[[i]] <- cbind(meth, cbind(gex, resp))
    }
    message(paste0("...->Plotting ",gene," (",site,"; ",region,") ","for ",drug_label," in ",ct))
    #print(to[[i]])
    #to <<- tos
    
    if(all(index_meth == "all")){
      p[[i]] <- to[[i]]
    }else{
      
      site <- sym(site)
      gene <- sym(gene)
      drug <- sym(drug)
      if(!use_tcga){
        target <- as.character(unique(drug_info$PUTATIVE_TARGET[drug_info$DRUG_NAME == drug_label]))
        if(gex_y){
          p[[i]] <- ggplot(data =to[[i]],
                           aes(x = !!site, y = !!(if(gex_y){gene}else{drug}))) + geom_point(aes(size=!!(if(gex_y){drug}else{gene})))+ # color = !!drug
            geom_smooth(method = "lm",color = "grey")+
            xlab(paste0("Methylation of ",site,"\nin the ",gene," gene ",region,""))+
            labs(y = if(gex_y){paste0("Gene expression\nof ",gene," in ",ct)}else{paste0("Drug response (AUC)\nto ",drug_label," (",target,")")},
                 size = if(gex_y){paste0("Drug response (AUC)\nto ",drug_label," (",target,")")}else{paste0("Gene expression\nof ",gene," in\n",ct)}
                 #shape = "Mutational status",
            )+ #color = paste0("Drug response (AUC)\nto ",drug_label,"")
            scale_size(trans = if(gex_y){'reverse'}else{"identity"})+
            #scale_color_viridis(option = "D")+
            theme_minimal()
        }else{
          p[[i]] <- ggplot(data =to[[i]],
                           aes(x = !!site, y = !!(if(gex_y){gene}else{drug}))) + geom_point(aes(color=!!(if(gex_y){drug}else{gene})))+ # color = !!drug
            geom_smooth(method = "lm",color = "grey")+
            xlab(paste0("Methylation of ",site,"\nin the ",gene," gene ",region,""))+
            labs(y = if(gex_y){paste0("Gene expression\nof ",gene," in ",ct)}else{paste0("Drug response (AUC)\nto ",drug_label," (",target,")")},
                 color = if(gex_y){paste0("Drug response (AUC)\nto ",drug_label," (",target,")")}else{paste0("Gene expression\nof ",gene," in\n",ct)}
                 #shape = "Mutational status",
            )+ #color = paste0("Drug response (AUC)\nto ",drug_label,"")
            scale_size(trans = if(gex_y){'reverse'}else{"identity"})+
            #scale_color_viridis(option = "D")+
            theme_minimal()+
            scale_color_viridis()
        }
      }else{
        p[[i]] <- ggplot(data =to[[i]],
                         aes(x = !!site, y = !!(if(gex_y){gene}else{gene}))) + geom_point()+ # color = !!drug
          geom_smooth(method = "lm",color = "grey")+
          xlab(paste0("Methylation of ",site,"\nin the ",gene," gene ",region,""))+
          labs(y = if(gex_y){paste0("Gene expression\nof ",gene," in ",ct)}else{paste0("Drug response (AUC)\nto ",drug_label,"")}
               #size = if(gex_y){paste0("Drug response (AUC)\nto ",drug_label,"")}else{paste0("Gene expression\nof ",gene," in ",ct)}
               #shape = "Mutational status",
          )+ #color = paste0("Drug response (AUC)\nto ",drug_label,"")
          scale_size(trans = if(gex_y){'reverse'}else{"identity"})+
          #scale_color_viridis(option = "D")+
          theme_minimal()
      }
    }
  }
  
  return(p)
}
################################################################################################


get_genomics <- function(RESULTS, # RESULTS-object
                         df = df, #all DMRS
                         which.cancer.type = which.cancer.type,
                         drivers = mutation
                         
){
  
  RESULTS$DAT$covariates <- RESULTS$DAT$covariates[row.names(RESULTS$DAT$bem),,drop = F]  #cbind(, RESULTS$DAT$bem[,drivers, drop = F])
  df_here <- df[df$cancertype == which.cancer.type,]
  df_here <- unnest(df_here[,c("Drug","cpg","index","gene","metht","methS","Gene")])
  df_here <- df_here[df_here$cpg %in% colnames(RESULTS$DAT$meth),]
  df_here <- df_here[df_here$methS >=3, ] #Filter, 3 is the same as used for filtering dDMRs
  # add gene expression correlation
  int <- intersect(row.names(RESULTS$DAT$meth),row.names(RESULTS$DAT_gex$meth))
  df_here$gex_cor <- lapply(1:nrow(df_here),
                            function(x) tryCatch(cor(RESULTS$DAT$meth[int,unlist(df_here[x,"cpg"])],RESULTS$DAT_gex$meth[int,unlist((df_here[x,"Gene"])[[1]][1])],use = "complete.obs"),error=function(e)NA)) %>% unlist
  df_here$gex_cor[is.na(df_here$gex_cor)] <- 0
  # make the before filter to have extended distributions
  df_here <- as.data.frame(df_here)
  df_here$atleastsome <- lapply(1:nrow(df_here), function(y){ # df_here$cpg
    screen <- RESULTS$DAT$resp #before::: #screen <- RESULTS$DAT$resp[, df_here[y,"Drug"],drop=F]
    screened_lines <- row.names(screen)[which(!is.na(screen))]
    screened_lines <- screened_lines[screened_lines %in% row.names(RESULTS$DAT$bem)]
    mutants <- names((RESULTS$DAT$bem[screened_lines,drivers])[RESULTS$DAT$bem[screened_lines,drivers]==1])
    wildtypes <- names((RESULTS$DAT$bem[screened_lines,drivers])[RESULTS$DAT$bem[screened_lines,drivers]==0])
    if(!df_here[y,"cpg"] %in% colnames(RESULTS$DAT$meth)){
      next
    }else{
      #under_m<-length(which(RESULTS$DAT$meth[mutants,df_here[y,"cpg"]]>0.8));
      #over_m<-length(which(RESULTS$DAT$meth[mutants,df_here[y,"cpg"]]<0.2));
      #under_w<-length(which(RESULTS$DAT$meth[wildtypes,df_here[y,"cpg"]]>0.8));
      #over_w<-length(which(RESULTS$DAT$meth[wildtypes,df_here[y,"cpg"]]<0.2));
      over <- length(which(RESULTS$DAT$meth[,df_here[y,"cpg"]]>0.7));
      under <- length(which(RESULTS$DAT$meth[,df_here[y,"cpg"]]<0.3));
    }
    #tmp <- min(na.omit(under_m),na.omit(over_m),na.omit(over_w),na.omit(under_w))
    tmp <- min(na.omit(under),na.omit(over))
    return(tmp)
  })%>%unlist
  
  df_here <- df_here[df_here$atleastsome>2,]# Final filter
  
  if(length(df_here$index%>%unique) == 0){cpg_asso_lrmut_full <- NULL;cpg_asso_lrmut <- NULL;}
  i <- 1; for(indices in df_here$index%>% unique){
    obj <- RESULTS$DAT
    obj$resp <- RESULTS$DAT$resp[,drivers,drop = F] #obj$resp <- RESULTS$DAT$resp[,df_here[df_here$index ==indices,"Drug"]%>%unique%>%unlist, drop=F]
    obj$meth <- RESULTS$DAT$meth[,df_here[df_here$index ==indices,"cpg"]%>%unlist, drop =F]
    print(colnames(obj$meth))
    skip <- ncol(obj$meth) # small hack
    #print(obj$meth)
    if(ncol(obj$meth==1)){obj$meth <- cbind(obj$meth,data.frame(none=rep(1, nrow(obj$meth))))}
    #obj$covariates <- cbind(obj$covariates, RESULTS$DAT$bem[,drivers, drop = F])
    obj <<- obj
    #cpg_asso_lrmut <- asso(data.object = obj, annotation.object = RESULTS$ANN, metadata = F, min_screened = 15, description = what, scaling = F, 
    #                       covariates = obj$covariates, min_auc = 0.7, min_responder = 3, lrtest = drivers)
    ##### ABOVE FUNCTION IS FOR DETECTING INTERACTIONS
    cpg_asso_lrmut <- asso_mut_meth(data.object = obj, annotation.object = RESULTS$ANN, metadata = F, min_screened = 15, description = what, scaling = F, 
                                    covariates = obj$covariates, min_auc = 1, min_responder = 1, lrtest = F)
    if(!is.null(cpg_asso_lrmut)){ # for dlbc, trametinib 1372-GDSC1 is in df, which was supposed to be filtered in pre-
      cpg_asso_lrmut$hit <- df_here[df_here$index ==indices,"gene"]%>%unique%>%unlist
      cpg_asso_lrmut$alt <- rep(drivers, nrow(cpg_asso_lrmut))
    }else{
      cpg_asso_lrmut <- data.frame(matrix(NA, nrow = skip, ncol = 11)); colnames(cpg_asso_lrmut) = c("P.Value","effectsize","t-value","FDR","CPG.Labels","indep","gene","amount","confidence","hit","alt")
    }
    if(i==1){cpg_asso_lrmut_full <- cpg_asso_lrmut}else{cpg_asso_lrmut_full <- rbind(cpg_asso_lrmut_full,cpg_asso_lrmut)}
    i <- i+1
  }; cpg_asso_lrmut_full <- cpg_asso_lrmut_full[cpg_asso_lrmut_full$gene != "none",]; #print(cpg_asso_lrmut_full)
  ## Exploiting the same order, please preserve, used instead of -log(sd) as confidence
  if(!is.null(cpg_asso_lrmut_full)){
    cpg_asso_lrmut_full$confidence <- abs(df_here$gex_cor)
    cpg_asso_lrmut_full$cor_gex <- df_here$gex_cor
    
    cpg_asso_lrmut_full$annotation <- cpg_asso_lrmut_full$alt
    cpg_asso_lrmut_full$gene <- cpg_asso_lrmut_full$hit
  }
  return(cpg_asso_lrmut_full)
}



### Function: Scatter plots for some hits
get_data_for_scatter <- function(df = df, results_collected = results_collected, which.cancer.type = which.cancer.type){ 
  df_fig2 <- df[df$cancertype == which.cancer.type,]
  RESULTS <- results_collected[[which.cancer.type]]
  RESULTS$DAT$gene_anno <- DAT$gene_anno # makes faster to just do the meth=T once
  RESULTS$DAT_gex$meth <- loadRData(paste0("results_deg/",
                                           which.cancer.type,".genes.RData"))$DAT_gex$meth
  RESULTS$gex_asso <- loadRData(paste0("results_deg/",
                                       which.cancer.type,".genes.RData"))$gex_asso
  
  dmrs_3 <- unnest(df_fig2, cols = c("cpg","enchancer","dhs","regulator","pvalues","effectsizes","genepos"))
  dmrs_3 <- dmrs_3[!is.na(dmrs_3$pvalues),]
  #dmrs_3$Drug <- unlist(lapply(dmrs_3$drug, function(x) ANN$drugs_info$DRUG_NAME[ANN$drugs_info$DRUG_ID == strsplit(x,"-")[[1]][1]][1]))
  dmrs_3$conn <- lapply(1:nrow(dmrs_3), function(i) paste0(dmrs_3$drug[i],"---",dmrs_3$Gene[i]))
  #dmrs_3$confidence <- unlist(lapply(dmrs_3$cpg, function(x){ ### dont need that here since it is done in the preparation of df
  #  if(!x %in% colnames(RESULTS$DAT$meth)){
  #    return(NA)
  #  }else{
  #    under<-length(which(RESULTS$DAT$meth[,x]>0.5));
  #    over<-length(which(RESULTS$DAT$meth[,x]<0.5));
  #  }
  #  return(min(under,over))
  #}))
  print(nrow(dmrs_3)) # same but the NA pvalues are removed
  return(list(dmrs_3=dmrs_3,RESULTS=RESULTS, which=which.cancer.type))
}
################################################################################################