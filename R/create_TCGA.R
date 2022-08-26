suppressPackageStartupMessages({
library(TCGAbiolinks)
library(dplyr)
library(DT)
library(ChAMP)
library(SummarizedExperiment)
})

create_query <- function(
### Creates a query of Illumina Human Methylation 450k and RNAseq
######################################################
  cancer.type = NULL, # what cancer type is used (in TCGA project)
  sample.type
){
  query.met <- GDCquery(project = paste("TCGA-", cancer.type, sep = ""),
                        data.category = "DNA Methylation",
                        legacy = FALSE,
                        platform = c("Illumina Human Methylation 450"),
                        sample.type = sample.type
                        )
  query.exp <- GDCquery(project = paste("TCGA-", cancer.type, sep = ""),
                        data.category = "Transcriptome Profiling",
                        data.type = "Gene Expression Quantification", 
                        workflow.type = "HTSeq - Counts",
                        sample.type = sample.type
                        )
  
  common.patients <- intersect(substr(getResults(query.met, cols = "cases"), 1, 12),
                               substr(getResults(query.exp, cols = "cases"), 1, 12)) # filter query for common patients
  
  query.met <- GDCquery(project = paste("TCGA-", cancer.type, sep = ""),
                        data.category = "DNA Methylation",
                        legacy = FALSE,
                        platform = c("Illumina Human Methylation 450"),
                        barcode = common.patients,
                        sample.type = sample.type
                        )
  query.exp <- GDCquery(project = paste("TCGA-", cancer.type, sep = ""),
                        data.category = "Transcriptome Profiling",
                        data.type = "Gene Expression Quantification", 
                        workflow.type = "HTSeq - Counts",
                        barcode = common.patients,
                        sample.type = sample.type
                        )
  
  return(list(exp = query.exp, met = query.met))
}######################################################


download_query <- function(
### Downloads the previously created query, adds a path to the object
######################################################
  query = NULL, # query object
  path.to.data = "TCGA/", # where the downloaded data is stored
  download.exp = TRUE, # should exp be downloaded or not
  download.met = TRUE # should met be downloaded or not
){
  if(unlist(query$exp$project) == unlist(query$met$project)){
    cancer.type <- unique(c(query$exp$project,query$met$project)) # define cancer type
  } else {
    print("The query projects do not match !")
  }
  path.to.data <- path.to.data
  dir.create(path.to.data, showWarnings = F)
  if(download.exp){
    print("Downloading expression query...")
    GDCdownload(query = query$exp, directory = path.to.data)
  }
  if(download.met){
    print("Downloading methylation query...")
    GDCdownload(query = query$met, directory = path.to.data, files.per.chunk = 10)
  }
  
  return(list(exp = query$exp, met = query$met, path = path.to.data))
}######################################################


create_data <- function(
### creates the data objects, meaning the matrices (creates raw betavalues)
######################################################
  query = NULL, # query object
  create.exp = TRUE, # should exp be made or not
  create.met = TRUE, # should met be made or not
  path.to.metadata = "metadata/TCGA/", # where the data objects should be stored
  load.metadata = FALSE # should metadata be loaded
){
  if(unlist(query$exp$project) == unlist(query$met$project)){
    cancer.type <- unique(c(query$exp$project,query$met$project)) # define cancer type
  } else {
    print("The query projects do not match !")
  }
  save(query, file = paste(path.to.metadata,unlist(cancer.type),"-data_q.RData",sep=""))
  
  path.to.data <- query$path
  if(create.exp & !load.metadata){
    print("Creating expression data...")
    data.exp <- GDCprepare(query$exp, directory = query$path)
    data.exp.path <- paste(path.to.metadata,unlist(cancer.type),"-data_exp_prepared_raw.RData",sep="")
    save(data.exp, file = data.exp.path)
  } else {
    if(create.exp){
      data.exp.path <- paste(path.to.metadata,unlist(cancer.type),"-data_exp_prepared_raw.RData",sep="")
      load(data.exp.path)
    } else {
      data.exp <- F
    }
  }
  if(create.met & !load.metadata){
    print("Creating methylation data...")
    data.met <- GDCprepare(query$met, directory = query$path)
    data.met.path <- paste(path.to.metadata,unlist(cancer.type),"-data_met_prepared_raw.RData",sep="")
    save(data.met, file = data.met.path)
  } else {
    if(create.met){
    data.met.path <- paste(path.to.metadata,unlist(cancer.type),"-data_met_prepared_raw.RData",sep="")
    load(data.met.path)
    } else {
      data.met <- F
    }
  }
  print("Data object is created, previous objects can be cleared !")
  return(list(exp = data.exp, met = data.met, cancer.type = unlist(cancer.type), path = query$path, query.exp = query$exp, query.met = query$met))
}######################################################


preprocess_DMR <- function(
### Downloads the previously created query
### Deprecated
######################################################
  data = NULL, # data object
  SampleCutoff = 0.5, # when to filter a sample (fraction of NA)
  ProbeCutoff = 0.05, # when to filter a probe (fraction of NA)
  method = "Delete", # imputation method
  phenotype = "shortLetterCode", # phenotype from data$met@colData$...
  path.to.metadata = "metadata/TCGA/" # where the data objects should be stored, does not exist yet
){
  data$exp <- NULL
  data$query.exp <- NULL # clear some space, because data object contains both methylation and GEx
  
  data.met.matrix <- SummarizedExperiment::assays(data$met)
  data$met@colData$Sample_Name <- data$met@colData$barcode
  myFilter <- champ.filter(beta = data.met.matrix@listData[[1]] %>% as.data.frame, pd = data$met@colData, arraytype = "450K", fixOutlier = F)
  data.met.matrix <- dim(data.met.matrix@listData[[1]]) # rm
  myImpute <- champ.impute(beta = myFilter$beta, pd = myFilter$pd, method = method, SampleCutoff = SampleCutoff, ProbeCutoff = ProbeCutoff)
  myFilter <- dim(myFilter) # rm
  myNorm <- champ.norm(beta = myImpute$beta, arraytype = "450K", cores = 1)
  pheno <- unlist(myImpute$pd[phenotype])
  if(phenotype == "shortLetterCode"){
    myNorm <- myNorm[, pheno %in% c("NT","TP")]
    pheno <- pheno[pheno %in% c("NT","TP")]
  }
  myImpute <- dim(myImpute) # rm
  gc()
  
  return(list(beta = myNorm, phenotype = pheno, phenotype_name = phenotype, cancer.type = data$cancer.type, query.met = data$met,
              summary = paste(
                paste("Raw Data: Dim ", as.character(data.met.matrix), sep = ""),
                paste("After filtering step: Dim ", myFilter, sep = ""),
                paste("After imputation step: Dim ", myImpute, sep = ""),
                paste("After normalization step: Dim ", dim(myNorm), sep = ""),
                sep = "\n")
              )
         )
}######################################################


calc_DEG <- function(
  ### Performs a differential expression analysis
  ### Deprecated
  ######################################################
  data = NULL, # data object
  path.to.metadata = "metadata/TCGA/", # where the data objects should be stored
  metadata = F, # if true, then metadata is loaded
  cancer.type = NULL # which cancer type if data is not supplied
){
  data$met <- NULL
  data$query.met <- NULL # clear some space, because data object contains both methylation and GEx
  if(!is.null(data)){
    cancer.type <- data$cancer.type
  }
  if(metadata == F){
    
    # Tumor vs Normal Samples
    samplesDown <- getResults(data$query.exp,cols = c("cases"))
    
    dataSmTP <- TCGAquery_SampleTypes(barcode = samplesDown,
                                      typesample = "TP")
    
    dataSmNT <- TCGAquery_SampleTypes(barcode = samplesDown,
                                      typesample = "NT")
    # Filter data
    queryDown <- GDCquery(project = data$cancer.type, 
                          data.category = "Transcriptome Profiling",
                          data.type = "Gene Expression Quantification", 
                          workflow.type = "HTSeq - Counts", 
                          barcode = c(dataSmTP, dataSmNT))
    
    #GDCdownload(query = queryDown)
    
    dataPrep1 <- GDCprepare(query = queryDown, directory = data$path)#, 
                            #save = TRUE, 
                            #save.filename = "TCGA_BRCA_HTSeq_Countds.rda")
    data.exp <- dataPrep1
    data.exp.path <- paste(path.to.metadata,unlist(cancer.type),"-data_exp.RData",sep="")
    save(data.exp, file = data.exp.path)
    
    dataPrep <- TCGAanalyze_Preprocessing(object = dataPrep1, 
                                          cor.cut = 0.6,
                                          datatype = "HTSeq - Counts")
    
    
    dataNorm <- TCGAanalyze_Normalization(tabDF = dataPrep,
                                          geneInfo = geneInfoHT,
                                          method = "gcContent")
    dataNorm <- TCGAanalyze_Normalization(tabDF = dataPrep,
                                          geneInfo = geneInfoHT,
                                          method = "geneLength")
    
    dataFilt <- TCGAanalyze_Filtering(tabDF = dataNorm,
                                      method = "quantile", 
                                      qnt.cut =  0.25)
  
    
    dataDEG <- TCGAanalyze_DEA(mat1 = dataFilt[,dataSmTP],
                                mat2 = dataFilt[,dataSmNT],
                                pipeline="limma",
                                Cond1type = "Normal",
                                Cond2type = "Tumor",
                                fdr.cut = 0.01 ,
                                logFC.cut = 1,
                                method = "glmLRT", ClinicalDF = data.frame())
    
    save(dataDEG, file = paste(path.to.metadata,cancer.type,"-data_met_DEG_","",".RData",sep=""))
  }else{
    load(paste(path.to.metadata,cancer.type,"-data_met_DEG_",method,".RData",sep=""))
  }
  print("Done exporting differential analysis !!!! :)")
  
  return(list(DEG = dataDEG,
              method = "limma",
              cancer.type = cancer.type
  )
  )
}######################################################


calc_DMR <- function(
  ### Deprecated
  ######################################################
  data = NULL, # data object
  path.to.metadata = "metadata/TCGA/", # where the data objects should be stored
  method = "bumphunter", # or cate,lasso,diffVar,
  metadata = F, # if true, then metadata is loaded
  cancer.type = NULL, # which cancer type if data is not supplied
  i = 100 # which drug should be tested
){
  if(!is.null(data) & !any(method %in% "bestdmr")){
    cancer.type <- data$cancer.type
  }
  
  if(metadata == F){
    if(method == "dmp"){
      dataDMR <- champ.DMP(beta = data$beta, pheno = data$phenotype, adjPVal = 0.25)
    }
    
    if(method == "bumphunter"){
      dataDMR <- champ.DMR(beta = data$beta, pheno = data$phenotype, method="Bumphunter", cores = 1)
    }
    if(method == "lasso"){
      dataDMR <- champ.DMR(beta = data$beta, pheno = data$phenotype, method="ProbeLasso", cores = 1) # champ.DMR nor champ.DMR_my seem to work for probelasso
    }
    if(method == "cate"){
      dataDMR <- champ.DMR_my(beta = data$beta, pheno = data$phenotype, method="DMRcate", cores = 1, cateMethod = "differential")
    }
    if(any(method %in% "diffVar")){
      dataDMR <- champ.DMR_my(beta = data$beta, pheno = data$phenotype, method="DMRcate", cores = 1, cateMethod = "diffVar")
    }
    
    if(any(method %in% "bestdmr")){
      cancer.type <- data$which
      dataDMR <- champ.DMR_my(beta = data$meth, pheno = data$resp_factor[,i], method="DMRcate", cores = 1, cateMethod = "differential")
    }
  
    save(dataDMR, file = paste(path.to.metadata,cancer.type,"-data_met_DMR_",method,".RData",sep=""))
  }else{
    load(paste(path.to.metadata,cancer.type,"-data_met_DMR_",method,".RData",sep=""))
  }
  print("Done exporting differential analysis !!!! :)")
  
  return(list(DMR = dataDMR,
              method = method,
              cancer.type = cancer.type
  )
  )
}######################################################
