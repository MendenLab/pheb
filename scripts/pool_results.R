#!/usr/bin/env Rscript
# patternclustering_gdsc.R

args <- as.numeric(commandArgs(trailingOnly = TRUE))
print(paste0("Running with arguments: ",""))
print(args)
source("R/script_params.R")
load(paste("data/cosmic.RData","",sep=''))
library(MultiAssayExperiment)
library(ELMER.data)
library(ELMER)
#---> generateing this here #load("metadata/summaries/df_patterns_updated.RData") # df$id, df_pattern, df_pattern_gex, df_pattern_hits, df_pattern_hits_validated

if(TRUE){ # make: df$id, df_pattern, df_pattern_gex, df_pattern_hits, df_pattern_hits_validated
       # Summarize elmer
  # get tcga patterns
  tcga_p <- list(); for(type in types){
    tcga_pattern <- tryCatch(get(load(file = paste0("metadata/elmer_onlytumors/",type,"/res.drug.regionTRUE.RData"))), error = function(e) NA)
    if(class(tcga_pattern)=="logical"){tcga_pattern <- NULL}
    tcga_p[[type]] <- rbind(do.call(rbind, lapply(tcga_pattern$distal, function(x) x[["pair_pos"]])),
                            do.call(rbind,lapply(tcga_pattern$distal, function(x) x[["pair_neg"]])),
                            do.call(rbind,lapply(tcga_pattern$promoter, function(x) x[["pair_pos"]])),
                            do.call(rbind,lapply(tcga_pattern$promoter, function(x) x[["pair_neg"]]))
    )
  }; tcga_p <- do.call(rbind, tcga_p)
  
  # get gdsc patterns
  gdsc_p <- list(); for(type in types){
    tcga_pattern <- tryCatch(get(load(file = paste0("metadata/elmer_cl//",type,"/res.drug.regionFALSE.RData"))), error = function(e) NA)
    if(class(tcga_pattern)=="logical"){tcga_pattern <- NULL}
    gdsc_p[[type]] <- rbind(do.call(rbind, lapply(tcga_pattern$distal, function(x) x[["pair_pos"]])),
                            do.call(rbind,lapply(tcga_pattern$distal, function(x) x[["pair_neg"]])),
                            do.call(rbind,lapply(tcga_pattern$promoter, function(x) x[["pair_pos"]])),
                            do.call(rbind,lapply(tcga_pattern$promoter, function(x) x[["pair_neg"]]))
    )
  }; gdsc_p <- do.call(rbind, gdsc_p)
  
  tcga_tf <- list(); for(type in types){
    tcga_pattern <- tryCatch(get(load(file = paste0("metadata/elmer_onlytumors/",type,"/res.drug.regionTRUE.RData"))), error = function(e) NA)
    if(class(tcga_pattern)=="logical"){tcga_pattern <- NULL}
    tcga_tf[[type]] <- rbind(do.call(rbind, lapply(tcga_pattern$distal, function(x) x[["pair_pos"]])),
                            do.call(rbind,lapply(tcga_pattern$distal, function(x) x[["pair_neg"]])),
                            do.call(rbind,lapply(tcga_pattern$promoter, function(x) x[["pair_pos"]])),
                            do.call(rbind,lapply(tcga_pattern$promoter, function(x) x[["pair_neg"]]))
    )
  }; tcga_tf <- do.call(rbind, tcga_tf)
  
  # get gdsc patterns
  gdsc_tf <- list(); for(type in types){
    tcga_pattern <- tryCatch(get(load(file = paste0("metadata/elmer_cl//",type,"/res.drug.regionFALSE.RData"))), error = function(e) NA)
    if(class(tcga_pattern)=="logical"){tcga_pattern <- NULL}
    gdsc_tf[[type]] <- rbind(do.call(rbind, lapply(tcga_pattern$distal, function(x) x[["pair_pos"]])),
                            do.call(rbind,lapply(tcga_pattern$distal, function(x) x[["pair_neg"]])),
                            do.call(rbind,lapply(tcga_pattern$promoter, function(x) x[["pair_pos"]])),
                            do.call(rbind,lapply(tcga_pattern$promoter, function(x) x[["pair_neg"]]))
    )
  }; gdsc_tf <- do.call(rbind, gdsc_tf)
  
  
  ### form df_pattern for the elmer associations
  df$id <- 1:nrow(df)
  df_pattern <- df
  tcga_eqtl <- list();for(i in 1:nrow(df_pattern)){
    tcga_eqtl[[i]] <- c(tcga_p$Symbol[which((tcga_p$Probe %in% df_pattern$cpg[i][[1]]) & tcga_p$cancertype == df_pattern$annotation[i] & tcga_p$drug == df_pattern$drug[i])])
  }
  tcga_eqtl_label <- unlist(lapply(tcga_eqtl, function(x) paste0(paste(x%>%unique, collapse = "/"),"")))
  gdsc_eqtl <- list();for(i in 1:nrow(df_pattern)){
    gdsc_eqtl[[i]] <- c(gdsc_p$Symbol[which((gdsc_p$Probe %in% df_pattern$cpg[i][[1]]) & gdsc_p$cancertype == df_pattern$annotation[i] & gdsc_p$drug == df_pattern$drug[i])]) # & gdsc_p$cancertype == df_pattern$annotation[i] & gdsc_p$drug == df_pattern$drug[i]
  }
  gdsc_eqtl_label <- unlist(lapply(gdsc_eqtl, function(x) paste0(paste(x%>%unique, collapse = "/"),"")))
  
  # get directionality
  gdsc_eqtl_dir <- list();for(i in 1:nrow(df_pattern)){
    gdsc_eqtl_dir[[i]] <- c(gdsc_p$direction[which((gdsc_p$Probe %in% df_pattern$cpg[i][[1]]) & gdsc_p$cancertype == df_pattern$annotation[i] & gdsc_p$drug == df_pattern$drug[i] )])
  }
  tcga_eqtl_dir <- list();for(i in 1:nrow(df_pattern)){
    tcga_eqtl_dir[[i]] <- c(tcga_p$direction[which((tcga_p$Probe %in% df_pattern$cpg[i][[1]]) & tcga_p$cancertype == df_pattern$annotation[i] & tcga_p$drug == df_pattern$drug[i])])
  }
  
  df_pattern$gdsc_eqtl_label <- gdsc_eqtl_label
  df_pattern$tcga_eqtl_label <- tcga_eqtl_label
  df_pattern$gdsc_eqtl <- gdsc_eqtl
  df_pattern$tcga_eqtl <- tcga_eqtl
  df_pattern$gdsc_eqtl_dir <- gdsc_eqtl_dir
  df_pattern$tcga_eqtl_dir <- tcga_eqtl_dir
  
  ## add functional region annotation
  distal.probes <- get.feature.probe(promoter = F, # or promoters=FALSE for distal enhancers
                                     genome = "hg19", 
                                     met.platform = "450K", 
                                     rm.chr = paste0("chr",c("X","Y")))%>%as.data.frame%>%row.names
  promoter.probes <- get.feature.probe(promoter = T, # or promoters=FALSE for distal enhancers
                                     genome = "hg19", 
                                     met.platform = "450K", 
                                     rm.chr = paste0("chr",c("X","Y")))%>%as.data.frame%>%row.names
  distal <- lapply(df_pattern$cpg, function(x) any(x %in% distal.probes))%>%unlist
  promoter <- lapply(df_pattern$cpg, function(x) any(x %in% promoter.probes))%>%unlist
  distal[distal] <- "distal"
  distal[distal != "distal"] <- ""
  promoter[promoter] <- "promoter"
  promoter[promoter !="promoter"] <- ""
  df_pattern$elmer_region <- unlist(lapply(1:length(promoter), function(x) paste0(promoter[x],"-",distal[x])))
  
  ### 182 hits that are translatable
  load("metadata/summaries/all_updated.RData") # <-import previous validation results
  df_pattern <- cbind(df_pattern,all[,c("id",
                    "discovery","validation")])
  
  ### fill in gene expression and correlations with the elmer analysed genes
  df_gex_red_elmer <- df_gex[(df_gex$indep %in% df$Drug)&(df_gex$gene %in% unique(unlist(df_pattern$gdsc_eqtl))),]
  #which.row.gex <- lapply(1:nrow(df_pattern), function(i){
  #  if(length(unlist(df_pattern[i,"gdsc_eqtl"]))>0){
  #                                      unlist(as.character(lapply(1:length(unlist(df_pattern[i,"gdsc_eqtl"])), function(j){
  #                                        which((df_gex_red_elmer$indep == df_pattern[i,"Drug"])&
  #                                       (df_gex_red_elmer$gene %in% unlist(df_pattern[i,"gdsc_eqtl"])[j])&
  #                                       df_gex_red_elmer$annotation == df_pattern[i,"cancertype"]
  #                                       )})))
  #  }else{""}
  # })
  
  which.row.gex <- loadRData(fileName = "metadata/sorted_gex_asso_elmer_onlytumors.RData")
  
  df_pattern$P.Value.gex.elmer <- 
    lapply(which.row.gex, function(y){df_gex_red_elmer[lapply(as.numeric(y), function(x) if(length(x)==""){NA}else{x})%>%unlist,"P.Value"]})
  df_pattern$effectsize.gex.elmer <- 
    lapply(which.row.gex, function(y){df_gex_red_elmer[lapply(as.numeric(y), function(x) if(length(x)==""){NA}else{x})%>%unlist,"effectsize"]})
  
  df_pattern_gex <- df_pattern
  df_pattern_gex$gdsc_eqtl_dir <- lapply(1:length(df_pattern_gex$gdsc_eqtl), function(x) df_pattern_gex$gdsc_eqtl_dir[[x]][!duplicated(df_pattern_gex$gdsc_eqtl[[x]])])
  df_pattern_gex$gdsc_eqtl <- lapply(df_pattern_gex$gdsc_eqtl, function(x) unique(x))
  df_pattern_gex$P.Value.gex.elmer <- lapply(df_pattern_gex$P.Value.gex.elmer, function(x) unique(x))
  df_pattern_gex$effectsize.gex.elmer <- lapply(df_pattern_gex$effectsize.gex.elmer, function(x) unique(x))
  library(tidyverse)
  df_pattern_gex <- df_pattern_gex[,!duplicated(colnames(df_pattern_gex))]
  df_pattern_gex <- unnest(df_pattern_gex, cols = c("id","eff","gdsc_eqtl","P.Value.gex.elmer","effectsize.gex.elmer","gdsc_eqtl_dir"))
  # 1443 gene drug - interactions
  
  # Hits
  df_pattern_hits <- df_pattern[df_pattern$gdsc_eqtl_label!="" & df_pattern$tcga_eqtl_label!="",]
  df_pattern_hits <- df_pattern_hits[lapply(1:nrow(df_pattern_hits), function(x) any(unlist(df_pattern_hits$gdsc_eqtl[x]) %in% unlist(df_pattern_hits$tcga_eqtl[x])))%>%unlist,]
  df_pattern_hits_validated <- df_pattern_hits[(0<(df_pattern_hits$discovery * df_pattern_hits$validation)) & !is.na(df_pattern_hits$discovery),]
    # get direction -> filter out not matching directions, only WDR46 here.. 55 left
  df_pattern_hits$tcga_genes <- lapply(1:nrow(df_pattern_hits), function(i) paste0("",unique(intersect(df_pattern_hits$gdsc_eqtl[i][[1]],df_pattern_hits$tcga_eqtl[i][[1]]))))
  com1 <- lapply(1:nrow(df_pattern_hits), function(i) paste0("",unique(df_pattern_hits$gdsc_eqtl_dir[i][[1]][df_pattern_hits$gdsc_eqtl[i][[1]] %in% df_pattern_hits$tcga_genes[i]])))
  com2 <- lapply(1:nrow(df_pattern_hits), function(i) paste0("",unique(df_pattern_hits$tcga_eqtl_dir[i][[1]][df_pattern_hits$tcga_eqtl[i][[1]] %in% df_pattern_hits$tcga_genes[i]])))
  df_pattern_hits <- df_pattern_hits[unlist(com1) == unlist(com2),]
  
  if(F){View(df_pattern_hits);View(df_pattern_hits_validated)}
  
  
  #######################################################################################################################################
  save(df, df_pattern, df_pattern_gex, df_pattern_hits, df_pattern_hits_validated, file = "metadata/summaries/df_patterns_updated_onlytumors.RData")
  #######################################################################################################################################
}






# list_cancertypes is gex data
# results_collected is meth data

###########
##
## This makes post analyses after building df with dDMR-search using scripts for DMR calling.
## Further data generated taht is imported here:
## elmer.R for meth-gex associations in cell lines and tumours
## CCLE_validation.R for validation results
## hits_genomics.R for genomics associations
###########
library(IlluminaHumanMethylation450kmanifest)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(tibble)
library(cowplot)
library(jsonlite)
library(purrr)
library(dplyr)
library(tidyverse)
library(data.table)
library(tidyr)
library(plyr)
library(ggplot2)
library(viridis)
library(grid)
library(gridExtra)
library(ggforce)
library(pals)
library(ComplexHeatmap)
load(paste("metadata/cosmic.RData","",sep=''))
types <- names(table(cosmic$tissue)[table(cosmic$tissue) > 15 & table(cosmic$tissue) <100])
load("metadata/summaries/df_patterns_updated_onlytumors.RData") # df$id, df_pattern, df_pattern_gex, df_pattern_hits, df_pattern_hits_validated
#cancertypes_meth <- cancertypes # rewrite because then cancertypes can be gene expression
source("R/script_params.R", local = TRUE)
script_params(vector = NULL, submission = F, args = 1) # for importing libraries

METADATA_TCGA <- "metadata/TCGA_onlytumors/"

if(TRUE){
  list_cancertypes <- list()
  for(i in 1:length(types)){
    which.cancer.type <- script_params(vector = types 
                                       , submission = T, args = i, parallel = F)
    if(F){
      ANN <- init_annotations(if(which.cancer.type == "COAD"){"COREAD"}else{which.cancer.type}, # TODO this function should also be able to take CCLE data
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
      list_cancertypes[[i]] <- DAT$meth
    }
    RESULTS <- loadRData(paste0("metadata/results_deg/",
                                which.cancer.type,".genes.RData"))
    list_cancertypes[[i]] <- RESULTS$DAT_gex$meth
    print(which.cancer.type)
  }
  
  
  types_cancertypes <- lapply(1:length(list_cancertypes), function(x)
  {rep(types[x],nrow(list_cancertypes[[x]])) }) %>% unlist
  for(i in 1:length(list_cancertypes)){
    if(i ==1)
    {cancertypes <- list_cancertypes[[1]]}
    else
    {cancertypes <- rbind(cancertypes, list_cancertypes[[i]])}
    print(i)
  }
}
################################################################################################



### BUILTS df and plots it,  df_gex below
### NEEDS: ANN, RESULTS
if(T){
  for(i in types){
    print(i)
    cpg_asso <- loadRData(paste0("metadata/results_dmp_new_combp/",i,".regions.RData"))
    RESULTS <- loadRData(paste0("metadata/results_deg/",
                                i,".genes.RData"))
    cpg_asso$effectsize <- unlist(lapply(cpg_asso$effectsizes, function(x) mean(na.omit(x))))
    cpg_asso$P.Value <- unlist(lapply(cpg_asso$pvalue, mean))
    cpg_asso$confidence <- -cpg_asso$nprobes
    cpg_asso$gene <- lapply(cpg_asso$gene, function(x) unique(x)[1])
    cpg_asso$Drug <- cpg_asso$drug
    cpg_asso$drug <- unlist(lapply(cpg_asso$drug, function(x) ANN$drugs_info$DRUG_NAME[ANN$drugs_info$DRUG_ID == strsplit(x,"-")[[1]][1]][1]))
    cpg_asso$annotation <- rep(i, nrow(cpg_asso))
    cpg_asso$metht <- lapply(1:length(cpg_asso$cpg), function(y){
      tmp <- c()
      drugresponse <- cpg_asso$Drug[y]
      for(i in 1:length(cpg_asso$cpg[y])){
        if(!(unlist(cpg_asso$cpg[y]))[i] %in% colnames(RESULTS$DAT$meth)){
          
        }else{
          under<-length(which(RESULTS$DAT$meth[!is.na(RESULTS$DAT$resp[,drugresponse]),(unlist(cpg_asso$cpg[y]))[i]]>0.7));
          over<-length(which(RESULTS$DAT$meth[!is.na(RESULTS$DAT$resp[,drugresponse]),(unlist(cpg_asso$cpg[y]))[i]]<0.3));
        }
        tmp[i] <- min(na.omit(under),na.omit(over))
      }
      return(tmp)
    })
    cpg_asso$methS <- unlist(lapply(cpg_asso$metht, median))
    #df <- cpg_asso[abs(cpg_asso$effectsize) <4,]
    
    if(i == "ALL"){
      df_new <- cpg_asso
    } else {
      df_new <- rbind(cpg_asso, df_new)
    }
  }
  ###!!!!!!!!!! df <- df_new
  
  # add cancer genes
  gene_list <- unlist(lapply(strsplit("ABL1 CCDC6 EIF1AX HIST1H2BD MED12 POLE SMARCB1 UPF3A ACO1 CCND1 EIF2S2 HIST1H3B MED23 POT1 SMC1A VHL ACVR1 CD1D ELF3 HIST1H4E MEN1 POU2AF1 SMC3 WASF3 ACVR1B CD58 EML4 HLA-A MET POU2F2 SMO WT1 ACVR2A CD70 EP300 HLA-B MGA PPM1D SMTNL2 XIRP2 ACVR2B CD79A EPAS1 HLA-C MLH1 PPP2R1A SNX25 XPO1 ADNP CD79B EPHA2 HNF1A MPL PPP6C SOCS1 ZBTB20 AJUBA CDC27 EPS8 HOXB3 MPO PRDM1 SOX17 ZBTB7B AKT1 CDC73 ERBB2 HRAS MSH2 PRKAR1A SOX9 ZFHX3 ALB CDH1 ERBB3 IDH1 MSH6 PSG4 SPEN ZFP36L1 ALK CDH10 ERCC2 IDH2 MTOR PSIP1 SPOP ZFP36L2 ALPK2 CDK12 ERG IKBKB MUC17 PTCH1 SPTAN1 ZFX AMER1 CDK4 ESR1 IKZF1 MUC6 PTEN SRC ZMYM3 APC CDKN1A ETNK1 IL6ST MXRA5 PTPN11 SRSF2 ZNF471 APOL2 CDKN1B EZH2 IL7R MYD88 PTPRB STAG2 ZNF620 ARHGAP35 CDKN2A FAM104A ING1 MYOCD QKI STAT3 ZNF750 ARHGAP5 CDKN2C FAM166A INTS12 MYOD1 RAC1 STAT5B ZNF800 ARID1A CEBPA FAM46C IPO7 NBPF1 RACGAP1 STK11 ZNRF3 ARID1B CHD4 FAT1 IRF4 NCOR1 RAD21 STK19 ZRSR2 ARID2 CHD8 FBXO11 ITGB7 NF1 RASA1 STX2 ARID5B CIB3 FBXW7 ITPKB NF2 RB1 SUFU ASXL1 CIC FGFR1 JAK1 NFE2L2 RBM10 TBC1D12 ATM CMTR2 FGFR2 JAK2 NIPBL RET TBL1XR1 ATP1A1 CNBD1 FGFR3 JAK3 NOTCH1 RHEB TBX3 ATP1B1 CNOT3 FLG KANSL1 NOTCH2 RHOA TCEB1 ATP2B3 COL2A1 FLT3 KCNJ5 NPM1 RHOB TCF12 ATRX COL5A1 FOSL2 KDM5C NRAS RIT1 TCF7L2 AXIN1 COL5A3 FOXA1 KDM6A NSD1 RNF43 TCP11L2 AXIN2 CREBBP FOXA2 KDR NT5C2 RPL10 TDRD10 AZGP1 CRLF2 FOXL2 KEAP1 NTN4 RPL22 TERT B2M CSDE1 FOXQ1 KEL NTRK3 RPL5 TET2 BAP1 CSF1R FRMD7 KIT NUP210L RPS15 TG BCLAF1 CSF3R FUBP1 KLF4 OMA1 RPS2 TGFBR2 BCOR CTCF GAGE12J KLF5 OR4A16 RPS6KA3 TGIF1 BHMT2 CTNNA1 GATA1 KLHL8 OR4N2 RREB1 TIMM17A BIRC3 CTNNB1 GATA2 KMT2A OR52N1 RUNX1 TNF BMPR2 CUL3 GATA3 KMT2B OTUD7A RXRA TNFAIP3 BRAF CUL4B GNA11 KMT2C PAPD5 SELP TNFRSF14 BRCA1 CUX1 GNA13 KMT2D PAX5 SETBP1 TOP2A BRCA2 CYLD GNAQ KRAS PBRM1 SETD2 TP53 BRD7 DAXX GNAS KRT5 PCBP1 SF3B1 TRAF3 C3orf70 DDX3X GNB1 LATS2 PDAP1 SGK1 TRAF7 CACNA1D DDX5 GNPTAB LCTL PDGFRA SH2B3 TRIM23 CALR DIAPH1 GPS2 LZTR1 PDSS2 SLC1A3 TSC1 CARD11 DICER1 GTF2I MAP2K1 PDYN SLC26A3 TSC2 CASP8 DIS3 GUSB MAP2K2 PHF6 SLC44A3 TSHR CBFB DNM2 H3F3A MAP2K4 PHOX2B SLC4A5 TTLL9 CBL DNMT3A H3F3B MAP2K7 PIK3CA SMAD2 TYRO3 CBLB EEF1A1 HIST1H1C MAP3K1 PIK3R1 SMAD4 U2AF1 CCDC120 EGFR HIST1H1E MAX PLCG1 SMARCA4 UBR5",' '),function(x) gsub(" ","",paste(x,' ',sep=''))))
  df <- df[df$methS >=3,]
  # add rest
  df$pathway <- unlist(lapply(as.character(df$drug), function(x) as.character(drugs_info$PATHWAY_NAME[as.character(drugs_info$DRUG_NAME) == x][1])))
  df$cancertype <- df$annotation
  #df <- df[abs(df$effectsize)<3,] why
  df$confidence <- df$nprobes
  df$Gene <- df$gene
  df$gene <- unlist(lapply(1:nrow(df), function(i) paste0(df$gene[i],"---",df$drug[i])))
  df$cancergene <- df$Gene %in% gene_list
  ncg <- read.csv("metadata/cancergenes_list.txt", sep = "\t")
  df$GGene <- unlist(lapply(df$Gene, function(x) x[1]))
  df$ncg <- df$GGene %in% unique(as.character(ncg$X711_Known_Cancer_Genes))
  #df$cor is added in hits_synergy since mining frames
  
  ## SAVE IMAGE save.image(file = "/storage/groups/cbm01/workspace/alexander.ohnmacht/BEST/.RData")
}
################################################################################################




### Expression dataframe
if(T){
  for(i in types){
    print(i)
    RESULTS <- loadRData(paste0("metadata/results_deg/",
                                i,".genes.RData"))
    RESULTS$gex_asso$annotation <- rep(i, nrow(RESULTS$gex_asso))
    if(i == "ALL"){
      df_gex <- RESULTS$gex_asso
    } else {
      df_gex <- rbind(RESULTS$gex_asso, df_gex)
    }
  }
}
################################################################################################


### Merge validation with df ###################################################################
res <- lapply(
  list.files(path = "metadata/methylDB", recursive = T, pattern = "res_updated_3", full.names = T),
  function(x){tmp<-as.data.frame(loadRData(x));tmp$cancertype<-rep(strsplit(x,"/")[[1]][10],nrow(tmp));return(tmp)})
ccle_data_summary <- do.call(rbind, res)
#ccle_data_summary$Gene <- lapply(make.names(row.names(ccle_data_summary)),function(x) strsplit(x,"\\.")[[1]][1])%>%unlist

df$id <- 1:nrow(df)
#mg$id <- mg$index_df
ccle_data_summary$id <- ccle_data_summary$index.in.df
all <- merge(df, ccle_data_summary, by=c("id"), all=T)
colnames(all) = make.unique(colnames(all)) #DEPRECATED#x:CCLE,y:GDSC,,1:TCGA, except for cor: cox.x:GDSC, cor.y:CCLE
all <- all[!is.na(all$id),]
ncg <- read.csv("metadata/cancergenes_list.txt", sep = "\t")
all$ncg <- all$GGene %in% unique(as.character(ncg$X711_Known_Cancer_Genes))
save(all, file = "metadata/summaries/all_updated.RData")
df_pattern$discovery <- all$discovery
df_pattern$validation <- all$validation # some cancer types were missing
df_pattern$validation_not_gdsc <- all$validation_not_gdsc
################################################################################################


### Make TCGA data frames
if(T){
  library(biomaRt)
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  genes = getLDS(attributes = c("hgnc_symbol"), values = colnames(tcga_gex_elmer) , mart = human, attributesL = c("ensembl_gene_id"), martL = human, uniqueRows=T)
  tcga_gex_collected <- list()
  tcga_meth_collected <- list()
  for(type in types){
    print(type)
    use_tcga <- T
    tcga_gex_elmer <- tryCatch(loadRData(paste(METADATA_TCGA, type,"_gex.RData",sep="")),error = function(e) NULL)
    tcga_meth_elmer <- tryCatch(loadRData(paste(METADATA_TCGA, type,"_betas.RData",sep="")),error = function(e) NULL)
    #if((!is.null(tcga_gex_elmer) & !is.null(tcga_meth_elmer) & ncol(tcga_gex_elmer)>0 & ncol(tcga_gex_elmer)>0) | !use_tcga){
    ### big if
    if(length(tcga_meth_elmer)<7){tcga_meth_elmer <- tryCatch(tcga_meth_elmer$beta%>%t, error = function(e) NULL)}
    if(use_tcga){
      names <- lapply(colnames(tcga_gex_elmer), function(x) paste0("", unique(genes$HGNC.symbol[genes$Gene.stable.ID == strsplit(x,"\\.")[[1]][1]]))[1])%>% unlist
      colnames(tcga_gex_elmer) <- names
      tcga_gex_elmer <- tcga_gex_elmer[,!duplicated(colnames(tcga_gex_elmer))]
      tcga_gex_elmer <- tcga_gex_elmer[,colnames(tcga_gex_elmer) %in% genes$HGNC.symbol]
    }
    #}
    tcga_gex_collected[[type]] <- tcga_gex_elmer
    tcga_meth_collected[[type]] <- tcga_meth_elmer
  }
  save(tcga_gex_collected,tcga_meth_collected, file = "metadata/summaries/tcga_data_symbols_onlytumors.RData")
}
################################################################################################








if(TRUE){ # SUMMARIZE elmer results

  
  ### look for distance to drug targets from hits ########################################################################################################################
  load("metadata/summaries/df_patterns_updated_onlytumors.RData") # df$id, df_pattern, df_pattern_gex, df_pattern_hits, df_pattern_hits_validated
  library(OmnipathR)
  
  interactions <- import_Omnipath_Interactions()
  
  OPI_g <- interaction_graph(interactions = interactions)
  #OPI_g <- string_db$get_graph()
  source("R/create_GDSC.R")
  source("scripts/paths.R", local = TRUE)
  ANN2 <- init_annotations(if(which.cancer.type == "COAD"){"COREAD"}else{which.cancer.type}, # TODO this function should also be able to take CCLE data
                           drug=T,
                           meth=F,
                           cfe=F, 
                           path.to.cosmic= COSMIC_PATH, 
                           path.to.drugs_original2 = GDSC2_SCREEN_PATH,# new data
                           path.to.drugs_original1= GDSC1_SCREEN_PATH, # old data
                           path.to.drugs_info = DRUG_INFO_PATH,
                           path.to.reference_original = METH_MANIFEST_PATH)
  df_pattern$target <- unlist(lapply(df$Drug, function(x) ANN2$drugs_info$PUTATIVE_TARGET[ANN2$drugs_info$DRUG_ID  == strsplit(x,"-")[[1]][1]]))
  
  # define data frame
  df_pattern_path <- df_pattern[df_pattern$id %in% df_pattern_hits$id,]
  
  #prepare for algorithm
  df_pattern_path$target <- lapply(df_pattern_path$target, function(x) paste("",strsplit(as.character(x),", ")[[1]], sep = ""))
  df_pattern_path$target[df_pattern_path$target == "Telomerase"] <- list(c("TERT", "TERC", "DKC1", "TEP1"))
  df_pattern_path$target[df_pattern_path$target == "Microtubule stabiliser"] <- list(c("MOAP1","MAP1", "MAP2", "MAP4", "tau","DCX"))
  df_pattern_path$target[df_pattern_path$target == "Antimetabolite (DNA & RNA)"] <- list(c("TYMS"))
  df_pattern_path$target[df_pattern_path$target == "dsDNA break induction"] <- list(c("LIG1", "LIG3"))
  df_pattern_path$target[df_pattern_path$target == "NAE"] <- list(c("NEDD8"))
  df_pattern_path$target[df_pattern_path$target == "HSP90"] <- list(c("HSP90AA1","HSP90AB1"))
  df_pattern_path$target[df_pattern_path$target %in% list(c("MEK1","MEK2"))] <- list(c("MAP2K1","MAP2K2"))
  df_pattern_path$target[df_pattern_path$target == "TAK1"] <- list(c("MAP3K7"))
  logical <- unlist(lapply(df_pattern_path$target, function(x) "PI3Kalpha" %in% unlist(x)))
  df_pattern_path$target[logical] <- lapply(df_pattern_path$target[logical], function(x) c("PIK3CA",unlist(x)))
  logical <- unlist(lapply(df_pattern_path$target, function(x) "PI3Kgamma" %in% unlist(x)))
  df_pattern_path$target[logical] <- lapply(df_pattern_path$target[logical], function(x) c("PIK3CG",unlist(x)))
  
  path_list <- list()
  for(i in 1:nrow(df_pattern_path)){
    p <- unlist(df_pattern_path$target[i])
    q <- unique(intersect(df_pattern_path$gdsc_eqtl[i][[1]],df_pattern_path$tcga_eqtl[i][[1]]))
    grid <- expand.grid(p,q)
    grid_names <- grid
    
    #For string
    #p <- string_db$map(grid, "Var1", removeUnmappedRows = F)$STRING_id; grid$Var1 <- p
    #q <- string_db$map(grid, "Var2", removeUnmappedRows = F)$STRING_id; grid$Var2 <- q
    ###########
    
    #printPath_es(shortest_paths(OPI_g,
    #                            from = as.character(grid[j,"Var2"]),
    #                            to = as.character(grid[j,"Var1"]),
    #                            output = 'epath'
    #                            )$epath[[1]],OPI_g)
    
    
    #SHORTEST PATHS ALGORITHM
    if(T){
      allshortestpaths <- list()
      source <- "Var1"
      target <- "Var2"
      for(j in 1:nrow(grid)){
        #allshortestpaths[[j]] <- tryCatch(all_shortest_paths(OPI_g,
        #                                     from = as.character(grid[j,"Var2"]),
        #                                     to = as.character(grid[j,"Var1"]))$res, error = function(e) list(NULL))
        allshortestpaths[[j]] <- tryCatch(k_shortest_yen(graph = OPI_g, 
                                                         src = as.character(grid[j,source]), 
                                                         dest = as.character(grid[j,target]), k = 10), error = function(e) list(NULL))
        for(k in 10:1){
          if(any(is.null(unlist(allshortestpaths[[j]])))){
            allshortestpaths[[j]] <- tryCatch(k_shortest_yen(graph = OPI_g, 
                                                             src = as.character(grid[j,source]), 
                                                             dest = as.character(grid[j,target]), k = k), error = function(e) list(NULL))
          }
        }
      }; allshortestpaths_save <- allshortestpaths
      
      if(length((unlist(allshortestpaths)))==0){
        source <- "Var2" #from biomarker to target: "Var2"
        target <- "Var1"
        for(j in 1:nrow(grid)){
          #allshortestpaths[[j]] <- tryCatch(all_shortest_paths(OPI_g,
          #                                     from = as.character(grid[j,"Var2"]),
          #                                     to = as.character(grid[j,"Var1"]))$res, error = function(e) list(NULL))
          allshortestpaths[[j]] <- tryCatch(k_shortest_yen(graph = OPI_g, 
                                                           src = as.character(grid[j,source]), 
                                                           dest = as.character(grid[j,target]), k = 10), error = function(e) list(NULL))
          for(k in 10:1){
            if(any(is.null(unlist(allshortestpaths[[j]])))){
              allshortestpaths[[j]] <- tryCatch(k_shortest_yen(graph = OPI_g, 
                                                               src = as.character(grid[j,source]), 
                                                               dest = as.character(grid[j,target]), k = k), error = function(e) list(NULL))
            }
          }
        }; allshortestpaths <- c(allshortestpaths,allshortestpaths_save)
      }
    }
    #########################
    
    
    # insert almost shortest paths here
    #allpaths <- all_simple_paths(OPI_g,
    #                             from = as.character(grid[j,"Var2"]),
    #                             to = as.character(grid[j,"Var1"]))$res
    #printPath_vs(allshortestpaths,OPI_g)
    
    grid_names
    OPI_here <- induced_subgraph(OPI_g, v = unique((allshortestpaths %>% unlist)))
    colors <- colnames(as_adjacency_matrix(OPI_here));
    colors[!colors %in% c(q,p)]<-"grey90"; 
    colors[colors %in% c(p)]<-"pink2";
    colors[colors %in% c(q)]<-"turquoise3";
    if(length(colors)!=0){
      plot(OPI_here, 
           edge.color = adjustcolor("grey80", alpha.f = 1),
           label.font = "arial",
           vertex.label.dist=2,
           vertex.label.color="black",
           vertex.color = adjustcolor(colors, alpha.f = 1)
      )
      path <- recordPlot()
      invisible(dev.off())
      path_list[[i]] <- path
    }else{
      path_list[[i]] <- ""
    }
    #Filter out first degree connections
    #interactions_filtered <- dplyr::filter(interactions, 
    #                                   !(source_genesymbol %in% unique(names(allshortestpaths %>% unlist))) ,
    #                                   !(target_genesymbol %in% unique(names(allshortestpaths %>% unlist))))
    #all_shortest_paths(interaction_graph(interactions = interactions_filtered),
    #                   from = as.character(grid[j,"Var2"]),
    #                   to = as.character(grid[j,"Var1"]))$res
    #i <- i+1
  }
  df_pattern_path$path <- path_list
  save(df_pattern_path, file = "metadata/summaries/df_patterns_updated_path_pretty_onlytumors.RData")
  
}
