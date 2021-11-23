#!/bin/bash

# Check if virtualenv exists
ENVS=$(conda env list | awk '{print $1}' )
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
if [[ $ENVS = *"pheb_r3.5.3"* ]]; then
   conda activate "pheb_r3.5.3"
   if [[ $CONDA_PREFIX = *"pheb_r3.5.3"* ]]; then
     echo "Successfully initialised pheb."
   fi;
else 
   echo "Virtual environment does not exist, installing ..."
   conda env create -f ${DIR}"/r3.5.3_env.yml"
   conda activate "pheb_r3.5.3"
   if [[ $CONDA_PREFIX = *"pheb_r3.5.3"* ]]; then
     echo "Successfully initialised pheb."
   fi;
fi;

# Check if R is correct and install CRAN packages
if [[ `command -v R` = *"pheb_r3.5.3"* ]]; then
   echo "R version is correct."
   
   # Install BiocManager
   Rscript -e 'install.packages("BiocManager", repos="https://cloud.r-project.org")'
   
   # Install remaining packages
   Rscript -e 'BiocManager::install(c("stringr","parallel","doParallel","foreach","tidyr","dplyr","tidyverse","DT","DESeq2","reshape2","batchtools","mixtools","getopt","optparse","RcppArmadillo","argparser","CpGassoc","ggplot2","ggrepel","readxl","ChAMP","SummarizedExperiment","biomaRt","CePa","methylKit","EnsDb.Hsapiens.v86","IlluminaHumanMethylation450kanno.ilmn12.hg19","GenomicRanges","ggpubr","RBGL","graph","","MultiAssayExperiment","ELMER","ELMER.data","lmtest","limma","minfi","IlluminaHumanMethylation450kmanifest","tibble","cowplot","jsonlite","rjson","purrr","data.table","plyr","viridis","grid","gridExtra","ggforce","pals","OmnipathR","Rcpp","glmnet","knitr","devtools"))'

   # Install CRAN packages
   ## Rscript -e 'install.packages("", repos="https://cloud.r-project.org")'

   # Install Cairo packages
   export LD_LIBRARY_PATH=$CONDA_PREFIX/lib
   export CAIRO_INCLUDE_PATH=$CONDA_PREFIX/include/cairo
   export CAIRO_CFLAGS=-I$CONDA_PREFIX/include/cairo
   export CAIRO_LIBS=$CONDA_PREFIX/lib/cairo
   conda install -c conda-forge r-cairo
   Rscript -e 'BiocManager::install(c("ComplexHeatmap")'

   # Install TCGAbiolinks
   Rscript -e 'Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true")'
   Rscript -e 'remotes::install_github(repo="https://github.com/BioinformaticsFMRP/TCGAbiolinks")'

   # Install comb-p
   conda install -yc bioconda combined-pvalues

fi;


