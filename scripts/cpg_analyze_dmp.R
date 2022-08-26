what <- "GDSC"

# Create annotation objects for GDSC drugs and methylation array
ANN <- init_annotations(if(which.cancer.type == "COAD"){"COREAD"}else{which.cancer.type},
                        drug=T,
                        meth=T,
                        cfe=F, #<-----------------------------------------------------------------------
                        path.to.cosmic= COSMIC_PATH, 
                        path.to.drugs_original2 = GDSC2_SCREEN_PATH,
                        path.to.drugs_original1= GDSC1_SCREEN_PATH,
                        path.to.drugs_info = DRUG_INFO_PATH,
                        path.to.reference_original = METH_MANIFEST_PATH)

# Create and preprocess dataframes of cell line methylation SINGLE PROBES and drug sensitivity data
DAT <- init_data(annotation.object = ANN,
                 path.to.response.folder = SCREEN_MY_PATH, 
                 path.to.observation.folder = GDSC_METH_PATH,
                 path.to.bem.folder = BEM_PATH, 
                 path.to.beta.folder = METH_BETA_PATH,
                 skip_anno = T,
                 cellline_info = "data/Cell_Lines_Detail_GDSC.xlsx",
                 easy = T) 

selected <- grepl(drugid, colnames(DAT$resp))
DAT$resp <- subset(DAT$resp, select = selected)

if(DMP_GDSC){
  which.drugs.index <- 1:ncol(DAT$resp)
  which.drugs.index <- na.omit(which.drugs.index[(apply(DAT$resp,2 , function(x) 
    length(which(!is.na(x))) > 15 &
      length(which(x < 0.7)) >= 3 )) 
  ]) # get the drugs with cell lines screened more than min_screened
  DAT$resp <- DAT$resp[,which.drugs.index]

  # Make association tests and annotate
  cpg_asso <- asso(data.object = DAT, annotation.object = ANN, metadata = F, min_screened = 15, description = what, scaling = F, 
                   covariates = DAT$covariates, min_auc = 0.7, min_responder = 3)
  # Export objects
  RESULTS$DAT$resp <- DAT$resp
  RESULTS$DAT$meth <- NULL
  RESULTS$cpg_asso <- cpg_asso
  RESULTS$DAT$gene_anno <- lapply(RESULTS$DAT$gene_anno, function(x) as.character(x))
  save(RESULTS, file = paste("metadata/results_dmp_new/",which.cancer.type,"/",
                             which.cancer.type,"_",as.character(RESULTS$method),"_meth_results_",as.character(1),".RData",sep=""))
  print("Calcs have all ended for DMRs!!!")
}


if(DEG_GDSC){ # for gene expression
  gex_original <- read.csv(file = PATH_GEX_GDSC,
                           header = TRUE, sep="\t")
  
  ### Preprocess
  gex <- gex_original[gex_original$GENE_SYMBOLS!="",]
  library("EnsDb.Hsapiens.v86")
  hsens <- EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86
  allGenes_gex <- ensembldb::select(hsens, keys = as.character(gex$GENE_SYMBOLS) , columns = c("ENTREZID", "SYMBOL"), keytype = "SYMBOL")
  
  GEX <- gex[ ,lapply(colnames(gex), function(x) head(strsplit(x,'\\.'))[[1]][2]) %in% RESULTS$ANN$cosmic$cosmic_id[RESULTS$ANN$cosmic$tissue == RESULTS$ANN$which]]
  rownames(GEX) = GEX[,1]
  GEX <- GEX[, c(-1,-2)]
  
  celllines <- unlist(lapply(colnames(GEX), function(x) strsplit(x,"\\.")[[1]][2])) %in% row.names(RESULTS$DAT$meth)
  #GEX <- GEX[DEG$DEG$gene,celllines] #<<<<<<<<
  colnames(GEX) = unlist(lapply(colnames(GEX), function(x) strsplit(x,"\\.")[[1]][2]))
  
  # Rewrite RESULTS
  RESULTS$DAT_gex <- RESULTS$DAT
  RESULTS$DAT_gex$meth <- GEX%>%t
  RESULTS$DAT_gex$gene_anno <- colnames(GEX%>%t)
  RESULTS$DMR <- NULL
  RESULTS$DAT_gex$method <- "gex"
  RESULTS$ANN$gene_label_generate <- F # inflating memory
  RESULTS$ANN$gene_list_generate <- F
  
  #save <- RESULTS$DAT_gex$resp
  RESULTS$DAT_gex$resp <- RESULTS$DAT_gex$resp
  gex_asso <- asso(data.object = RESULTS$DAT_gex, annotation.object = RESULTS$ANN, metadata = F, min_screened = 15, description = "gex", scaling = F, logit.transform = F,
                   covariates = RESULTS$DAT_gex$covariates, min_auc = 0.7, min_responder = 3)
  # Export objects
  RESULTS$gex_asso <- gex_asso
  
  save(RESULTS, file = paste("metadata/results_deg/",which.cancer.type,"_",RESULTS$method,"_gex_results_",as.character(args[1]),".RData",sep=""))
  print("Calcs have all ended for DEGs!!!")
}

DAT_gex=NULL
gex_asso <- NULL
DMR <- NULL
cpg_asso <- NULL
RESULTS <- list(DMR=DMR, ANN=ANN, DAT=DAT, cpg_asso=cpg_asso, gex_asso=gex_asso, DAT_gex= DAT_gex, method= what)
RESULTS$ANN$gene_label_generate <- F
RESULTS$ANN$gene_list_generate <- F


