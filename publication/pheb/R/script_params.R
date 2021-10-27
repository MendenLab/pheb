script_params <- function(
### Sets the initial paramters of the running analysis script
######################################################
  vector, # which cancer types
  submission, # T/F if calcs locally or cluster
  args,
  parallel = F
){ 
  source("R/generic_functions.R")
  source("R/create_TCGA.R")
  source("R/create_GDSC.R")
  source("R/calc_asso.R")
  champ.DMR_my <<- dget("R/champ.DMR_my.R")
  cpg.assoc_my <<- dget("R/cpg.assoc_my.R")
  submission <<- submission
  suppressPackageStartupMessages({
  library(doParallel)
  library(foreach)
  library(tidyr)
  library(DESeq2)
  })
  
  # for new submission:
  args <- args[1]
  cores <- detectCores()
  if(parallel){
    cl <<- makeForkCluster(cores[1]-3)
  } else {
    cl <<- makeCluster(1)
  }
  registerDoParallel(cl)
  
  if(submission){ # Collection of cancer types to be investigated
    print(args)
    if(length(args)==0){
      print("No arguments supplied.")
    }else{
      project <- vector[args]
    }
    which.cancer.type <<- project
    source("metadata/paths.R")
  } else {
    which.cancer.type <- vector[args]
  }
  return(which.cancer.type)
}######################################################



