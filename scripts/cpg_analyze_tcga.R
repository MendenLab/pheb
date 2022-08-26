TCGA <- T
METADATA_TCGA <- "metadata/TCGA_onlytumors/"

if(TCGA){
  # Create query from GDC server
  object <- create_query(cancer.type = if(which.cancer.type == "COREAD"){"COAD"}else{which.cancer.type},
                         sample.type = c("Primary Tumor","Metastatic")
                         )
  # Download data from GDC server
  TCGA_PATH <- paste("metadata/TCGA/",as.character(which.cancer.type),"/",sep="")
  object <- download_query(query = object, 
                           path.to.data = TCGA_PATH, 
                           download.exp = F, 
                           download.met = F)
  data.object <- create_data(query = object,
                             create.exp = T,
                             create.met = T,
                             path.to.metadata = METADATA_TCGA, 
                             load.metadata = T)
# set load.metadata=T and create...= F, such that raw files are loaded from the TCGA metadata
  
  if(T){
    print("Process gex data...")
    # save gex data
    dds <- DESeqDataSet(data.object$exp, design = ~ 1)
    keep <- rowSums(counts(dds)) >= 10
    dds <- dds[keep,]
    
    dds <- vst(dds) # rlog takes too long
    gex_tcga <- t(assay(dds)) %>% as.data.frame
    save(gex_tcga, file = paste(METADATA_TCGA, which.cancer.type,"_gex.RData",sep=""))
    #######################################
    rm(object)
  }
    
  if(T){ 
    print("Process meth data...")
    # save meth data
    data.object.processed <- preprocess_DMR(data = data.object, path.to.metadata = METADATA_TCGA)
    if(!file.exists(paste(METADATA_TCGA, which.cancer.type,"_betas.RData",sep=""))){
      save(data.object.processed, file = paste(METADATA_TCGA, which.cancer.type,"_betas.RData",sep=""))
    } # put that on top
    #rm(data.object)
    data <- data.object.processed
  }
}

