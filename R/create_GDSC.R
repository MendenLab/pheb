### libraries
suppressPackageStartupMessages({
library(readxl)
})

init_annotations <- function(
### Initialization of Annotation Object (containing manifests and functions) for METHYLATION
######################################################
  what.cancer.type, 
  path.to.cosmic="metadata/cosmic.RData", 
  path.to.reference_original="data/GDSC_methylation/manifests/HumanMethylation450_15017482_v1-2.csv",
  path.to.drugs_original2="data/v17.3_fitted_dose_response.xlsx",
  path.to.drugs_original1="data/GDSC_drugs/PANCANCER_IC_Thu Aug  1 13_41_49 2019.csv", #add path to second cohort of drugs
  path.to.drugs_info="data/Drug_listSun Jun 23 18_27_32 2019.csv",
  drug=T,
  meth=T,
  cfe=F
){

  which <- what.cancer.type
  print(paste("### Annotation of cancer type",which,sep=" "))
  load(paste(path.to.cosmic,"",sep=''))
  cancer_types <- as.data.frame(sort(table(cosmic$tissue), decreasing = TRUE))
  
  if(meth){
    print("Load CpG annotations...")
    reference_original <- read.csv(path.to.reference_original, sep=',', header = TRUE, skip = 7)    #<------ fishy, preprocessing the annotation file here
    reference_features <- reference_original[-(485578:nrow(reference_original)),] # discard control
    reference_features <- reference_features[c("IlmnID", "UCSC_RefGene_Name", "UCSC_CpG_Islands_Name")]
    reference_features <- reference_features[1:485512,] # discard the SNPs whatever
    
    ### Annotation tables ####
    gene_label_generate <- function(island_names){ # make the label for annotation table for gene-cpg association
      tmp <- lapply(island_names,function(x) reference_features$UCSC_RefGene_Name[reference_features$UCSC_CpG_Islands_Name == as.character(x)])
      tmp <- lapply(tmp, function(y) paste(unlist(lapply(y, function(x) strsplit(as.character(x), ";"))),""))
      return(lapply(tmp, function(x) paste('g: ',toString(as.character(unique(x))))))
    }
    gene_list_generate <- function(island_names){ # make the sublist of annotation table for gene-cpg association
      tmp <- lapply(island_names,function(x) reference_features$UCSC_RefGene_Name[reference_features$UCSC_CpG_Islands_Name == as.character(x)])
      tmp <- lapply(tmp, function(y) paste(unlist(lapply(y, function(x) strsplit(as.character(x), ";"))),""))
      return(lapply(tmp, function(x) as.character(unique(x))))
    }
    ##########################
  }
  
  if(drug){                                                                                                             #<====== needs annotation of MY_ID
    print("Load drug annotations...")
    #drugs_original <- as.data.frame(read_excel(path.to.drugs_original1))
    drugs_original <- "TODO: MUST BE ADDED"
    drugs_info <- read.csv(path.to.drugs_info)
  }
  
  if(cfe){
    print("Load mutations annotations...")
    ### Mutation Annotation File
  }
  
  
  return(list(which=which,
              cosmic=cosmic,
              gene_label_generate=if(meth){gene_label_generate}else{FALSE},
              gene_list_generate=if(meth){gene_list_generate}else{FALSE},
              drugs_original=if(drug){drugs_original}else{FALSE},
              drugs_info=if(drug){drugs_info}else{FALSE},
              reference_features=if(meth){reference_features}else{FALSE},
              cfe=if(cfe){TRUE}else{FALSE})
  )
}######################################################


init_data <- function(
### Initialization of Molecular and Response Data
######################################################  
  path.to.response.folder="data/response_DRUG_ID_",
  path.to.observation.folder='data/observation_island_',
  path.to.bem.folder='data/CellLines_Mo_BEMs/',
  path.to.beta.folder='data/GDSC_methylation/',
  annotation.object,
  skip_anno = F,
  method = "Delete", # imputation method
  SampleCutoff = 0.5, # when to filter a sample (fraction of NA)
  ProbeCutoff = 0.1, # when to filter a probe (fraction of NA)
  cellline_info = "data/Cell_Lines_Detail_GDSC.xlsx",
  easy = T
){
  which <- annotation.object$which
  ### Import databases and load metadata
  print("Importing datasets...")
  response <- loadRData(paste(path.to.response.folder,"",which,".RData",sep=""))
  response <- response
  
  if(file.exists(paste(path.to.bem.folder,which,'_simple_MOBEM.rdata.tsv',sep=''))){
    bem_original <- read.table(paste(path.to.bem.folder,which,'_simple_MOBEM.rdata.tsv',sep=''), sep="\t", header=T, check.names = F)
  } else { # if bem does not exist, assign an half-empty matrix
    bem_original <- matrix(0, ncol = 1, nrow = nrow(response))
    row.names(bem_original) = row.names(response)
  }
  bem <- as.data.frame(bem_original)
  
  cellline_info <- read_excel("data/Cell_Lines_Detail_GDSC.xlsx")%>%as.data.frame
  ### Preprocess
  print("Preprocess...")
  if(!easy){ ### deprecated

  }else{
    
    betas <- loadRData(paste(path.to.observation.folder,"methyl_processed_",which,".RData",sep=""))
    
  }
  
  meth <- as.data.frame(betas) # make as dataframe
  resp <- as.data.frame(response)
  
  if(file.exists(paste(path.to.bem.folder,which,'_simple_MOBEM.rdata.tsv',sep=''))){ # if bem exists
    row.names(bem) = bem[,1]
    bem <- t(bem[,-1]) # discard the rowname-column
    bem <- bem[, order(colSums(bem),decreasing = T)] # sort cfe for high abundancy
    cfe_type <- sapply(strsplit(sapply(strsplit(as.character(colnames(bem)), "_"),tail,1),":"),head,1) # type of cfe 
    #bem <- bem[,cfe_type %in% c("gain","loss","mut")] # filter out methylation
    bem <- bem[,cfe_type != "HypMET"]
  }
  
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
  
  ### Make covariates
  covariates <- cellline_info[which(cellline_info$`COSMIC identifier` %in% row.names(betas)),]
  row.names(covariates) = covariates$`COSMIC identifier`
  covariates <- covariates[,c("Growth Properties","Screen Medium","Microsatellite \r\ninstability Status (MSI)")]
  for(i in 1:ncol(covariates)){
    covariates[,i] <- as.factor(covariates[,i])
  }
  covariates <- covariates[,as.logical(apply(covariates,2, function(x) length(levels(factor(x)))) >1),drop=F]
  
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
              bem=if(any(annotation.object$cfe !=F)){bem}else{FALSE},
              covariates=covariates)
  )
}######################################################
