#!/usr/bin/env Rscript


args <- as.numeric(commandArgs(trailingOnly = TRUE))
print(paste0("Running with arguments: ",""))
print(args)
source("R/script_params.R")
load(paste("metadata/cosmic.RData","",sep=''))
use_tcga <- TRUE # either perform ELMER analysis for GDSC or TCGA

library(MultiAssayExperiment)
library(ELMER.data)
library(ELMER)
library(dplyr)
library(tidyr)
library(tidyverse)
library(biomaRt)

### Customize types
types <- names(table(cosmic$tissue)[table(cosmic$tissue) > 15 & table(cosmic$tissue) <100])
which.cancer.type <- script_params(vector = types, # cancer types
                                   submission = T, 
                                   args = args, 
                                   parallel = F) 
###################


type <- types[as.numeric(args)] # plus how many already done
message(paste0("STARTED: ",type,"(number:", type,")"))
print(paste0("STARTED: ",type,"(number:", type,")"))

human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
genes = getLDS(attributes = c("hgnc_symbol"), values = colnames(tcga_gex_elmer) , mart = human, attributesL = c("ensembl_gene_id"), martL = human, uniqueRows=T)
# ALL,KIRC,MM have no methylation data
# BLCA has no columns in gene expression data
# LAML has no rows in methylation data
# NB,SCLC have no methylation data -> not in TCGA

if(T){
  which.cancer.type <- type
  #which.cancer.type <- "COREAD"
  S <- get_data_for_scatter(df = df, results_collected = results_collected, 
                            which.cancer.type = which.cancer.type)
 basedir <- paste0("metadata/elmer",if(use_tcga){"_onlytumors"}else{"_cl"},"/",which.cancer.type,"/")
 dir.create(basedir, showWarnings = T)
  #Import TCGA data
  METADATA_TCGA <- "metadata/TCGA_onlytumors/" #"metadata/TCGA/"
  tcga_gex_elmer <- tryCatch(loadRData(paste(METADATA_TCGA, which.cancer.type,"_gex.RData",sep="")),error = function(e) NULL)
  tcga_meth_elmer <- tryCatch(loadRData(paste(METADATA_TCGA, which.cancer.type,"_betas.RData",sep="")),error = function(e) NULL)
  if((!is.null(tcga_gex_elmer) & !is.null(tcga_meth_elmer) & ncol(tcga_gex_elmer)>0 & ncol(tcga_gex_elmer)>0) | !use_tcga){
    ### big if
    if(use_tcga){
      if(length(tcga_meth_elmer)<7){tcga_meth_elmer <- tcga_meth_elmer$beta%>%t}
      
      # not needed for some reason, look at it
      #names <- lapply(colnames(tcga_gex_elmer), function(x) paste0("", unique(genes$Gene.stable.ID[genes$HGNC.symbol == strsplit(x,"\\.")[[1]][1]]))[1])%>% unlist
      #colnames(tcga_gex_elmer) <- names
      tcga_gex_elmer <- tcga_gex_elmer[,!duplicated(colnames(tcga_gex_elmer))]
      
      tcga_gex_elmer <- tcga_gex_elmer[,colnames(tcga_gex_elmer) %in% genes$Gene.stable.ID]
      
    }
  
    if(!use_tcga){ # GDSC meth and gex manual data
      met <- S$RESULTS$DAT$meth%>%t
      exp <- S$RESULTS$DAT_gex$meth%>%t
      
      # use mart
      names <- lapply(row.names(exp), function(x) paste0("", unique(genes$Gene.stable.ID[genes$HGNC.symbol == strsplit(x,"\\.")[[1]][1]]))[1])%>% unlist
      row.names(exp) <- names
      exp <- exp[!duplicated(row.names(exp)),intersect(row.names(S$RESULTS$DAT$meth),row.names(S$RESULTS$DAT_gex$meth))]
      exp <- exp[row.names(exp) %in% genes$Gene.stable.ID,intersect(row.names(S$RESULTS$DAT$meth),row.names(S$RESULTS$DAT_gex$meth))]
      met <- met[,intersect(row.names(S$RESULTS$DAT$meth),row.names(S$RESULTS$DAT_gex$meth))]
      met <- met[,intersect(row.names(S$RESULTS$DAT$meth),row.names(S$RESULTS$DAT_gex$meth))]
      colnames(exp) <- make.names(colnames(exp))
      colnames(met) <- make.names(colnames(met))
      
      assay <- c(rep("DNA methylation", ncol(met)),
                 rep("Gene expression", ncol(exp)))
      primary <- c(colnames(met),colnames(exp))
      colname <- c(colnames(met),colnames(exp))
      sampleMap <- data.frame(assay,primary,colname)
      
      # get covariates
      cov <- S$RESULTS$DAT$covariates
      row.names(cov) = make.names(row.names(cov))
      cov <- cov[colnames(met),,drop = F]
      #
      
      colData <- cov
      colData$primary <- row.names(colData)
      colData$type <- "tumour"
      colData$type <- factor(colData$type)
    }
    
    # Define promoters
    res.drug.region <- list()
    for(region in c("promoter","distal")){
      print(paste0("...Run for: ","(region: ", region,")"))
      distal.probes <- get.feature.probe(promoter = if(region=="promoter"){T}else{F}, # or promoters=FALSE for distal enhancers
                                         genome = "hg19", 
                                         met.platform = "450K", 
                                         rm.chr = paste0("chr",c("X","Y")))
      
      # Make data object
      if(use_tcga){ # for TCGA data
        mae <- createMAE(exp = if(use_tcga){tcga_gex_elmer%>%t}else{exp}, 
                       met = if(use_tcga){tcga_meth_elmer%>%t}else{met},
                       save = F,
                       linearize.exp = F,
                       save.filename = "mae.rda",
                       filter.probes = distal.probes,
                       met.platform = "450K",
                       genome = "hg19",
                       #colData = if(use_tcga){NULL}else{colData},
                       #sampleMap = if(use_tcga){NULL}else{sampleMap},
                       TCGA = use_tcga)
      }else{ # try manual data
        mae <- createMAE(exp = if(use_tcga){tcga_gex_elmer%>%t}else{exp}, 
                         met = if(use_tcga){tcga_meth_elmer%>%t}else{met},
                         save = F,
                         linearize.exp = F,
                         save.filename = "mae.rda",
                         filter.probes = distal.probes,
                         met.platform = "450K",
                         genome = "hg19",
                         colData = if(use_tcga){NULL}else{colData},
                         sampleMap = if(use_tcga){NULL}else{sampleMap},
                         TCGA = use_tcga)
      }
      
      #stratify for drug
      drugsall <- levels(factor(S$dmrs_3$drug))
      res.drug <- list()
      for(drug in drugsall){
        print(paste0("   ...Run for: ","(drug: ", drug,")"))
        # Define DDMRs
        used_probes <- S$dmrs_3$cpg[S$dmrs_3$drug == drug]
        #used_probes_granges <- makeGRangesFromDataFrame(S$dmrs_3[S$dmrs_3$cpg %in% used_probes & S$dmrs_3$drug == "Trametinib",]%>%column_to_rownames("cpg"), keep.extra.columns = T)
        
        nearGenes <- tryCatch(GetNearGenes(data = mae, 
                                  probes = used_probes, 
                                  numFlankingGenes = 20), # 10 upstream and 10 dowstream genes
                              error = function(e) NULL)
        
        pair_neg <- tryCatch(get.pair(data = mae, 
                              group.col = if(!use_tcga){"type"}else{"definition"},
                              group1 =  if(!use_tcga){"tumour"}else{"Primary solid Tumor"},
                              group2 = if(!use_tcga){""}else{""}, #Solid Tissue Normal
                              nearGenes = nearGenes,
                              mode = "unsupervised", 
                              permu.dir = paste0(basedir,"permu","-",region,"-neg","---",drug),
                              correlation = "negative",
                              permu.size = 50000, # Please set to 100000 to get significant results
                              raw.pvalue = 0.05,   
                              Pe = 0.001, # Please set to 0.001 to get significant results
                              filter.probes = TRUE, # See preAssociationProbeFiltering function
                              filter.percentage = 0.05,
                              filter.portion = 0.3,
                              dir.out = paste0(basedir,"pairs","-neg"),
                              cores = 5,
                              label = paste0(region,"---",drug),
                              save = T),
                          error = function(e) NULL)
        pair_pos <- tryCatch(get.pair(data = mae, 
                                  group.col = if(!use_tcga){"type"}else{"definition"},
                                  group1 =  if(!use_tcga){"tumour"}else{"Primary solid Tumor"},
                                  group2 = if(!use_tcga){""}else{""}, #Solid Tissue Normal
                                  nearGenes = nearGenes,
                                  mode = "unsupervised",
                                  permu.dir = paste0(basedir,"permu","-",region,"-pos","---",drug),
                                  correlation = "positive",
                                  permu.size = 50000, # Please set to 100000 to get significant results
                                  raw.pvalue = 0.05,   
                                  Pe = 0.001, # Please set to 0.001 to get significant results
                                  filter.probes = TRUE, # See preAssociationProbeFiltering function
                                  filter.percentage = 0.05,
                                  filter.portion = 0.3,
                                  dir.out = paste0(basedir,"pairs","-pos"),
                                  cores = 5,
                                  label = paste0(region,"---",drug),
                                  save = T),
                         error = function(e) NULL)
        
        #simulation does not work for non-tcga data :()
        if(!is.null(pair_neg)){if(nrow(pair_neg)>0){
          pair_neg$drug <- drug
          pair_neg$region <- region
          pair_neg$cancertype <- which.cancer.type
          pair_neg$direction <- "negative"
        }}
        
        if(!is.null(pair_pos)){if(nrow(pair_pos)>0){
          pair_pos$drug <- drug
          pair_pos$region <- region
          pair_pos$cancertype <- which.cancer.type
          pair_pos$direction <- "positive"
        }}
        
        enriched.motif_neg <- tryCatch(get.enriched.motif(data = mae,
                                             probes = pair_neg$Probe, 
                                             dir.out = paste0(basedir,"emotifs","-neg"), 
                                             label = paste0(region,"--neg--",drug),
                                             min.incidence = 3,
                                             lower.OR = 1.1), error = function(e) NULL)
        
        enriched.motif_pos <- tryCatch(get.enriched.motif(data = mae,
                                                      probes = pair_pos$Probe, 
                                                      dir.out = paste0(basedir,"emotifs","-pos"), # potentiall add ,"-pos") in the end of calcs
                                                      label = paste0(region,"--pos--",drug),
                                                      min.incidence = 3,
                                                      lower.OR = 1.1), error = function(e) NULL)
        
        TF_neg <- tryCatch(get.TFs(data = mae, 
                      group.col = if(!use_tcga){"type"}else{"definition"},
                      group1 =  if(!use_tcga){"tumour"}else{"Primary solid Tumor"},
                      group2 = if(!use_tcga){""}else{""},
                      mode = "unsupervised",
                      enriched.motif = enriched.motif_neg,
                      dir.out = paste0(basedir,"tfs"), 
                      cores = 1, 
                      label = paste0(region,"--neg--",drug)), error = function(e) NULL)
        TF_pos <- tryCatch(get.TFs(data = mae, 
                                            group.col = if(!use_tcga){"type"}else{"definition"},
                                            group1 =  if(!use_tcga){"tumour"}else{"Primary solid Tumor"},
                                            group2 = if(!use_tcga){""}else{""},
                                            mode = "unsupervised",
                                            enriched.motif = enriched.motif_pos,
                                            dir.out = paste0(basedir,"tfs"), 
                                            cores = 1, 
                                            label = paste0(region,"--pos--",drug)), error = function(e) NULL)
        
        # Accumulate data for drugs
        res.drug[[drug]] <- list(
          pair_pos = pair_pos,
          pair_neg = pair_neg,
          enriched.motif_pos = enriched.motif_pos,
          enriched.motif_neg = enriched.motif_neg,
          TF_pos = TF_pos,
          TF_neg = TF_neg
        )
      } # end if both meth and gex exist
      
      
      res.drug.region[[region]] <- res.drug
    }
    # Accumulate data for regions
      
  }else{res.drug.region <- NA}
  save(res.drug.region, file = paste0(basedir,"res.drug.region",as.character(use_tcga),".RData"))
  # Accumulate data for cancer types
  print(paste0("FINISHEDCT: ",which.cancer.type,"(number:", which.cancer.type,")"))
  message(paste0("FINISHEDCT: ",which.cancer.type,"(number:", which.cancer.type,")"))
} #cancer type end



