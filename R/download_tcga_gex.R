library(biomaRt)
for(which.cancer.type in (types[1:22])){
  if(file.exists(paste(METADATA_TCGA, which.cancer.type,"_gex.RData",sep=""))){
        gex_tcga <- loadRData(paste(METADATA_TCGA, which.cancer.type,"_gex.RData",sep=""))
        
        human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
        genes = getLDS(attributes = c("ensembl_gene_id"), values = colnames(gex_tcga) , mart = human, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
        colnames(gex_tcga) <- lapply(colnames(gex_tcga), function(x) paste0("", unique(genes$HGNC.symbol[genes$Gene.stable.ID == strsplit(x,"\\.")[[1]][1]]))[1])%>% unlist
        gex_tcga <- gex_tcga[,colnames(gex_tcga)!="" & !(duplicated(colnames(gex_tcga)))]
        
        
      }else{
        # Create query from GDC server
        object <- tryCatch(create_query(cancer.type = if(which.cancer.type == "COREAD"){"COAD"}else{which.cancer.type}), error = function(e) "noproject")
        # Download data from GDC server
        if(object != "noproject"){
        object <- download_query(query = object, 
                               path.to.data = paste0("TCGA/",if(which.cancer.type == "COREAD"){"COAD"}else{which.cancer.type}), 
                               download.exp = T, 
                               download.met = F) #<-----------------------------------------------------------------------
        data.object <- create_data(query = object, 
                                 create.exp = T, 
                                 create.met = F, 
                                 path.to.metadata = METADATA_TCGA, 
                                 load.metadata = F) #<----------------------------------------------------------
      
        dds <- DESeqDataSet(data.object$exp, design = ~ 1)
        keep <- rowSums(counts(dds)) >= 10
        dds <- dds[keep,]
        dds <- vst(dds) # rlog takes too long
        gex_tcga <- t(assay(dds)) %>% as.data.frame
        }
      }
    save(gex_tcga, file = paste(METADATA_TCGA, which.cancer.type,"_gex.RData",sep=""))
}

