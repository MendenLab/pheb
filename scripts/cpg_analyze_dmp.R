#!/usr/bin/env Rscript

args <- as.numeric(commandArgs(trailingOnly = TRUE))
print(paste0("Running with arguments: ",""))
print("Arguments:")
print(args)


source("R/script_params.R")
load(paste("metadata/cosmic.RData","",sep=''))
types <- names(table(cosmic$tissue)[table(cosmic$tissue) > 15 & table(cosmic$tissue) <100])
which.cancer.type <- script_params(vector = types, # cancer types
                                   submission = T, 
                                   args = args, 
                                   parallel = F) 


print(which.cancer.type)
load(paste(COSMIC_PATH, "",sep=''))
DMP_GDSC <- T
fill <- F # specifies if inital run is done and checks and calculates missing
numberofargs <- 10
what <- "GDSC"

METADATA_TCGA <- "metadata/TCGA_onlytumors/"


# Create annotation objects for GDSC drugs and methylation array
ANN <- init_annotations(if(which.cancer.type == "COAD"){"COREAD"}else{which.cancer.type}, # TODO this function should also be able to take CCLE data
                        drug=T,
                        meth=T,
                        cfe=F, #<-----------------------------------------------------------------------
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



DAT_gex=NULL
gex_asso <- NULL
DMR <- NULL
cpg_asso <- NULL
RESULTS <- list(DMR=DMR, ANN=ANN, DAT=DAT, cpg_asso=cpg_asso, gex_asso=gex_asso, DAT_gex= DAT_gex, method= what)
RESULTS$ANN$gene_label_generate <- F # inflating memory
RESULTS$ANN$gene_list_generate <- F

if(DEG_GDSC){
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
  RESULTS$DAT_gex$resp <- cut_df(RESULTS$DAT_gex$resp, numberofargs, args)
  gex_asso <- asso(data.object = RESULTS$DAT_gex, annotation.object = RESULTS$ANN, metadata = F, min_screened = 15, description = "gex", scaling = F, logit.transform = F,
                   covariates = RESULTS$DAT_gex$covariates, min_auc = 0.7, min_responder = 3)
  # Export objects
  RESULTS$gex_asso <- gex_asso
  
  save(RESULTS, file = paste("metadata/results_deg/",which.cancer.type,"_",RESULTS$method,"_gex_results_",as.character(args[1]),".RData",sep=""))
  print("Calcs have all ended for DEGs!")
}