#!/usr/bin/env Rscript


source("R/script_params.R")
library(dplyr)
library(readxl)
library(lmtest)
load(paste("metadata/cosmic.RData","",sep=''))
load("metadata/summaries/df_patterns_updated_path_pretty_onlytumors.RData") # df$id, df_pattern, df_pattern_gex, df_pattern_hits, df_pattern_hits_validated
dim(df)
types <- names(table(cosmic$tissue)[table(cosmic$tissue) > 15 & table(cosmic$tissue) <100])
script_params(vector = types, submission = T, args = 1, parallel = F)


if(TRUE){ # needs to only run once 
  results_collected <- list()
  for(type in types){
        print(type)
    ANN <- init_annotations(if(which.cancer.type == "COAD"){"COREAD"}else{type}, # TODO this function should also be able to take CCLE data
                            drug=T,
                            meth=F,
                            cfe=T,
                            path.to.cosmic= COSMIC_PATH, 
                            path.to.drugs_original2 = GDSC2_SCREEN_PATH,# new data
                            path.to.drugs_original1= GDSC1_SCREEN_PATH, # old data
                            path.to.drugs_info = DRUG_INFO_PATH,
                            path.to.reference_original = METH_MANIFEST_PATH)
    
    # Create and preprocess dataframes of cell line methylation SINGLE PROBES and drug sensitivity data
    DAT <- init_data(annotation.object = ANN,
                     path.to.response.folder = SCREEN_MY_PATH, #<----updated drug cohort, files are in pharmaco4
                     path.to.observation.folder = GDSC_METH_PATH, #change file that GDSC in front
                     path.to.bem.folder = BEM_PATH, 
                     path.to.beta.folder = METH_BETA_PATH,
                     skip_anno = T,
                     cellline_info = "data/Cell_Lines_Detail_GDSC.xlsx",
                     easy = T)
    
    for(i in 1:15){
        results_collected[[type]] <- tryCatch(loadRData(
          paste("metadata/results_dmp/","/","/",
        type,"_","GDSC","_meth_results_",as.character(i),".RData",sep="")), error = function(e) NULL)
        if(!is.null(results_collected[[type]])){break}
    }
    results_collected[[type]]$DAT$bem <- DAT$bem
  }
  ANN <- init_annotations(if(which.cancer.type == "COAD"){"COREAD"}else{type}, # TODO this function should also be able to take CCLE data
                          drug=T,
                          meth=T,
                          cfe=F,
                          path.to.cosmic= COSMIC_PATH, 
                          path.to.drugs_original2 = GDSC2_SCREEN_PATH,# new data
                          path.to.drugs_original1= GDSC1_SCREEN_PATH, # old data
                          path.to.drugs_info = DRUG_INFO_PATH,
                          path.to.reference_original = METH_MANIFEST_PATH)
  
  # Create and preprocess dataframes of cell line methylation SINGLE PROBES and drug sensitivity data
  DAT <- init_data(annotation.object = ANN,
                   path.to.response.folder = SCREEN_MY_PATH, #<----updated drug cohort, files are in pharmaco4
                   path.to.observation.folder = GDSC_METH_PATH, #change file that GDSC in front
                   path.to.bem.folder = BEM_PATH, 
                   path.to.beta.folder = METH_BETA_PATH,
                   skip_anno = T,
                   cellline_info = "data/Cell_Lines_Detail_GDSC.xlsx",
                   easy = T)
}



##### MAKE BEM TESTS
int_gen <- list()
for(type in types){
  which.cancer.type <- type
  print(which.cancer.type)
  ### MAKE lrtest
#which.cancer.type <- "LIHC"
  RESULTS <- results_collected[[which.cancer.type]]
  RESULTS$DAT_gex <- loadRData(paste0("metadata/results_deg/",
                                                              which.cancer.type,".genes.RData"))$DAT_gex
  RESULTS$DAT$gene_anno <- DAT$gene_anno # makes faster to just do the meth=T once
  
  mutations <- tryCatch(names(apply(RESULTS$DAT$bem, 2, sum)[apply(RESULTS$DAT$bem, 2, sum)>5]), error = function(e) NULL)
  if(!is.null(mutations)){
    gresults <- list()
    
    RESULTS$DAT$resp <- RESULTS$DAT$bem
    
    for(i in 1:length(mutations)){
      gresults[[i]] <- get_genomics(RESULTS = RESULTS, 
                                        df = df, 
                                        which.cancer.type = which.cancer.type, 
                                        drivers = mutations[i])
    };gresults <- do.call("rbind", gresults)
    
    if(!is.null(gresults)){
      gresults$fdr <- p.adjust(gresults$P.Value, method = "BH")
      #gresults <- gresults[abs(gresults$effectsize)<3,]
      gresults$cancertype <- rep(which.cancer.type, nrow(gresults))
      }
    int_gen[[which.cancer.type]] <- gresults
    save(gresults,file = paste0("metadata/results_genomics/asso_list-",which.cancer.type,"_v2.RData"))
  }else{
    int_gen[[which.cancer.type]] <- NULL
  }
}
save(int_gen,file = "metadata/results_genomics/asso_list_mut_meth_v2.RData")


if(TRUE){
  
  gresults <- do.call(rbind, loadRData("metadata/results_genomics/asso_list_mut_meth.RData"))
  gresults <- gresults[!is.na(gresults$P.Value),]
  gresults$annotation <- gresults$cancertype
  
  gresults$identifier <- unlist(lapply(1:nrow(gresults), function(i) paste0(gresults$gene[i],"---",gresults$alt[i])))
  gresults <- aggregate(list(P.Value = gresults$P.Value, effectsize = gresults$effectsize,  confidence = gresults$confidence, fdr = gresults$fdr), by = list(hit = factor(gresults$identifier), gene = gresults$identifier, id = gresults$gene, annotation = gresults$annotation, cancertype = gresults$cancertype, indep = gresults_old$indep, alt = gresults_old$alt), mean )
  gresults$cpgs <- lapply(gresults$hit, function(x) gresults_old$CPG.Labels[gresults_old$identifier == x])
  gresults$fdr_overall <- p.adjust(gresults$P.Value, method = "BH")
  gresults_old <- gresults
  gresults <- gresults[gresults$id %in% df_pattern_hits$gene,] #FILTER
  gresults$fdr_overall <- p.adjust(gresults$P.Value, method = "BH")
  gresults$Drug <- unlist(lapply(gresults$id, function(x) (df$Drug[df$gene %in% x])[1]))
  gresults$GGene <- unlist(lapply(gresults$id, function(x) (df$GGene[df$gene %in% x])[1]))
  p1 <- plot_volcano(df = gresults, 
               sign = gresults$fdr_overall < 0.1, 
               Cond = gresults$fdr_overall < 0.1, 
               annotation.object = RESULTS$ANN,
               annotations = "")
  plot_volcano_2(df = gresults,max_eff = 1.5, max_p = 11, eff = 0.2, pval = 0.1,
                 Cond = gresults$fdr_overall >= 0.1, 
                 annotation.object = RESULTS$ANN,
                 annotations = "")+
    labs(color = "Cancer type")+
    scale_size_continuous(name = "correlation with\ngene expression (abs)")+
    ggtitle("")
  

  

}






