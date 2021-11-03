init_annotations_preprocess <- function(
  ### Initialization of Annotation Object (containing manifests and functions) for METHYLATION
  ######################################################
  path.to.cosmic="metadata/cosmic.RData", 
  #path.to.reference_original="../pharmaco1/data_original/GDSC_methylation/manifests/HumanMethylation450_15017482_v1-2.csv", # use GDSC meth
  path.to.drugs_original="data/CCLE/",
  drug=T,
  meth=T,
  cfe=F
){
  temp <- list()
  load(paste(path.to.cosmic,"",sep=''))
  cancer_types <- as.data.frame(sort(table(cosmic$tissue), decreasing = TRUE))

  
  # Load drugs_info and create the metadata files for response
  drugs_info_ccle <- read.csv(paste(path.to.drugs_original,"CTRPv2.0_2015_ctd2_ExpandedDataset/v20.meta.per_compound.txt",sep=""), sep = "\t")
  drugs_info_ccle$DRUG_ID <- drugs_info_ccle$master_cpd_id
  drugs_info_ccle$DRUG_NAME <- drugs_info_ccle$cpd_name
  cl_info_ccle <- read.csv(paste(path.to.drugs_original,"CTRPv2.0_2015_ctd2_ExpandedDataset/v20.meta.per_cell_line.txt",sep=""), sep = "\t")
  experiment_ccle <- read.csv(paste(path.to.drugs_original,"CTRPv2.0_2015_ctd2_ExpandedDataset/v20.meta.per_experiment.txt",sep=""), sep = "\t")
  cl_anno_ccle <- read.csv(paste(path.to.drugs_original,"Cell_lines_annotations_20181226.txt",sep=""), sep = "\t")
  sample_info <- read.csv(paste0(path.to.drugs_original,"response/sample_info.csv"))
  
  if(drug){
    print("Make drug dataframes...")
    ccle_responses <- read.csv(paste(path.to.drugs_original,"CTRPv2.0_2015_ctd2_ExpandedDataset/","v20.data.curves_post_qc.txt",sep=""), sep = "\t")
    ccle_responses <- ccle_responses[,c("experiment_id","master_cpd_id","area_under_curve")]
    ccle_responses <- as.data.frame(pivot_wider(data = ccle_responses, names_from = master_cpd_id, values_from = area_under_curve))
    cl_master_id <- make.unique(as.character(unlist(lapply(ccle_responses$experiment_id, function(x) unique(experiment_ccle$master_ccl_id[experiment_ccle$experiment_id == x]) ))))
    row.names(ccle_responses) = cl_master_id
    ccle_responses <- ccle_responses[,-1]
    cl_info_ccle$cosmic_id <- as.character(lapply(cl_info_ccle$ccl_name, function(x)  sample_info$COSMIC_ID[x == sample_info$stripped_cell_line_name]))
    cl_info_ccle$cosmic_id[cl_info_ccle$cosmic_id == "NA"] <- NA
    cl_info_ccle$cosmic_id[cl_info_ccle$cosmic_id == "numeric(0)"] <- NA
    cl_info_ccle$cosmic_id <- as.character(lapply(cl_info_ccle$ccl_name, function(x)  sample_info$COSMIC_ID[x == sample_info$stripped_cell_line_name]))
    cl_info_ccle$tissue_gdsc <- as.character(lapply(cl_info_ccle$cosmic_id, function(x)  na.omit(unique(cosmic$tissue[x == cosmic$cosmic_id & !is.na(x)]))))
    cl_info_ccle$tissue <- as.character(lapply(cl_info_ccle$cosmic_id, function(x)  na.omit(unique(cosmic$tissue[x == cosmic$cosmic_id & !is.na(x)]))))
    for(name in names(table(cl_info_ccle$tissue_gdsc))){ # matching histology and primary site with not classified cell lines
      tmp1 <- as.character(unique(cl_info_ccle$ccle_primary_site[cl_info_ccle$tissue_gdsc == name]))
      tmp1 <- tmp1[tmp1 != ""]
      tmp2 <- as.character(unique(cl_info_ccle$ccle_primary_hist[cl_info_ccle$tissue_gdsc == name]))
      #tmp2 <- tmp2[tmp2 != ""]
      tmp3 <- as.character(unique(cl_info_ccle$ccle_hist_subtype_1[cl_info_ccle$tissue_gdsc == name]))
      #tmp3 <- tmp3[tmp3 != ""]
      
      if(length(tmp1) == 1){
        cl_info_ccle$tissue[cl_info_ccle$ccle_primary_site == tmp1 & (cl_info_ccle$ccle_primary_hist %in% tmp2 & cl_info_ccle$ccle_hist_subtype_1 %in% tmp3)] <- name
      }
      print("Make drug response table for:")
      print(name)
      celllines <- unlist(lapply(cl_info_ccle$master_ccl_id[cl_info_ccle$tissue == name], function(x) experiment_ccle$experiment_id[experiment_ccle$master_ccl_id == x]))
      celllines_overlap <- unlist(lapply(cl_info_ccle$master_ccl_id[cl_info_ccle$tissue_gdsc == name], function(x) experiment_ccle$experiment_id[experiment_ccle$master_ccl_id == x]))
      temp[[name]] <- ccle_responses[celllines,]
      #save(temp, file = paste("data/CCLE/response_CCLE_",name,".RData",sep=""))
    }
  }
  
  if(cfe){
    print("Load mutations annotations...")
    ### Mutation Annotation File
  }
  
  
  return(list(
              cosmic=cosmic,
              resp=temp,
              gene_label_generate=if(meth){}else{FALSE},
              gene_list_generate=if(meth){}else{FALSE},
              drugs_original=if(drug){}else{FALSE},
              drugs_info=drugs_info_ccle,
              cl_info= cl_info_ccle,
              reference_features=if(meth){}else{FALSE},
              cfe=if(cfe){TRUE}else{FALSE})
  )
}######################################################





init_data_preprocess <- function(
  ### Initialization of Molecular and Response Data
  ######################################################  
  path.to.response.folder="data/CCLE/response_CCLE_",
  path.to.observation.folder='observation_island_',
  path.to.bem.folder='data/CellLines_Mo_BEMs/',
  path.to.beta.folder='data/CCLE/meth/',
  annotation.object,
  skip_anno = F,
  method = "Delete", # imputation method
  SampleCutoff = 0.5, # when to filter a sample (fraction of NA)
  ProbeCutoff = 0.1, # when to filter a probe (fraction of NA)
  easy = F # here the methylation gets live processed
){
  which <- annotation.object$which
  ### Import databases and load metadata
  print("Importing datasets...")
  response <- loadRData(paste(path.to.response.folder,"",which,".RData",sep=""))
  response <- response
  
  #bem_original <- read.table(paste(path.to.bem.folder,which,'_simple_MOBEM.rdata.tsv',sep=''), sep="\t", header=T, check.names = F)
  #bem <- as.data.frame(bem_original)
  
  ### Preprocess
  print("Preprocess...")
  if(!easy){
    
    ### initiate beta value matrix
    ccle_absolute_combined <- read_excel(paste(path.to.beta.folder,"CCLE_ABSOLUTE_combined_20181227.xlsx",sep=""))
    ccle_cgi <- read.csv(paste(path.to.beta.folder,"CCLE_RRBS_cgi_CpG_clusters_20181119.txt",sep=""), sep = "\t")
    ccle_tss1kb <- read.csv(paste(path.to.beta.folder,"CCLE_RRBS_TSS1kb_20181022.txt",sep=""), sep = "\t")
    ccle_enh <- read.csv(paste(path.to.beta.folder,"CCLE_RRBS_enh_CpG_clusters_20181119.txt",sep=""), sep = "\t")
    ccle_tss <- read.csv(paste(path.to.beta.folder,"CCLE_RRBS_tss_CpG_clusters_20181022.txt",sep=""), sep = "\t")
    
  }else{
    
    betas <- NULL
    
  }
  
  meth <- as.data.frame(betas) # make as dataframe
  resp <- as.data.frame(response)
  
  row.names(bem) = bem[,1]
  bem <- t(bem[,-1]) # discard the rowname-column
  bem <- bem[, order(colSums(bem),decreasing = T)] # sort cfe for high abundancy
  cfe_type <- sapply(strsplit(sapply(strsplit(as.character(colnames(bem)), "_"),tail,1),":"),head,1) # type of cfe 
  #bem <- bem[,cfe_type %in% c("gain","loss","mut")] # filter out methylation
  bem <- bem[,cfe_type != "HypMET"]
  
  if(any(annotation.object$drugs_info !=F) & any(annotation.object$reference_features !=F)){
    meth <- meth[intersect(rownames(resp),rownames(meth)),]
    resp <- resp[intersect(rownames(resp),rownames(meth)),]
  }
  
  if(any(annotation.object$cfe !=F) & any(annotation.object$reference_features !=F)){
    bem <- bem[intersect(rownames(bem),rownames(meth)),] 
    meth <- meth[intersect(rownames(bem),rownames(meth)),]
  }
  
  if(any(annotation.object$cfe !=F) & any(annotation.object$drugs_info !=F)){
    bem <- bem[intersect(rownames(bem),rownames(resp)),] 
    resp <- resp[intersect(rownames(bem),rownames(resp)),]
  }
  
  ### Annotation Table
  if(skip_anno==F){
    if(any(annotation.object$reference_features !=F)){
      print("Make Gene Annotations for CpG...")
      gene_label <- annotation.object$gene_label_generate(colnames(meth))
      gene_list <- annotation.object$gene_list_generate(colnames(meth))
    }
  }else{gene_label <- F; gene_list <- F}
  
  return(list(which=which,
              gene_label=if(any(annotation.object$reference_features !=F)){gene_label}else{FALSE},
              gene_list=if(any(annotation.object$reference_features !=F)){gene_list}else{FALSE},
              meth=if(any(annotation.object$reference_features !=F)){meth}else{FALSE},
              resp=if(any(annotation.object$drugs_info !=F)){resp}else{FALSE},
              bem=if(any(annotation.object$cfe !=F)){bem}else{FALSE})
  )
}######################################################