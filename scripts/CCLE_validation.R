#!/usr/bin/env Rscript

args <- as.numeric(commandArgs(trailingOnly = TRUE))
print(paste0("Running with arguments: ",""))
print(args)
source("R/script_params.R")
load(paste("metadata/cosmic.RData","",sep=''))


# CCLE validation
path <- "data/CCLE/"
source("R/create_CCLE.R")
source("R/script_params.R")
library(tidyverse)
chr_names <- read.csv("metadata/chr_names.txt", sep = "\t", header = F) %>% # for mapping NCBI to UCSC chromosomes
  mutate(V1 = as.character(V1)) %>% mutate(V2 = as.character(V2)); chr_names$V2[chr_names$V2 == ""] <- "none"
dim(df)

### Prepare gene expression of CCLE ###
if(!file.exists("data/CCLE/CCLE_RNAseq_processed.RData")){
  library(DESeq2)
  library(CePa)
  library(biomaRt)
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  gex_ccle <- read.gct("data/CCLE/CCLE_RNAseq_genes_counts_20180929.gct")
  genes = getLDS(attributes = c("ensembl_gene_id"), values = row.names(gex_ccle) , mart = human, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
  row.names(gex_ccle) <- lapply(row.names(gex_ccle), function(x) paste0("", unique(genes$HGNC.symbol[genes$Gene.stable.ID == strsplit(x,"\\.")[[1]][1]]))[1])%>% unlist
  gex_ccle <- gex_ccle[row.names(gex_ccle)!="" & !(duplicated(row.names(gex_ccle))),]
  names <- lapply(colnames(gex_ccle), function(x) strsplit(x,"_")[[1]][1]) %>% unlist
  colnames(gex_ccle) <- lapply(names, function(x) paste0("",ccle$cl_info$meth[ccle$cl_info$ccl_name == x]))%>%unlist
  gex_ccle <- gex_ccle[,colnames(gex_ccle)!=""] %>% t %>% as.data.frame %>% t
  ct_gex_ccle <- unlist(lapply(colnames(gex_ccle), function(x) paste0(ccle$cl_info$tissue[ccle$cl_info$meth ==x])))
  dds <- DESeqDataSetFromMatrix(countData = gex_ccle,
                                colData = as.data.frame(ct_gex_ccle),
                                design = ~ 1)
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep,]
  dds <- vst(dds) # rlog takes too long
  gex_ccle <- t(assay(dds)) %>% as.data.frame
  gex_ccle$cancertype <- ct_gex_ccle
  test <- gex_ccle %>% rownames_to_column() %>% group_split(cancertype)
  names <- unlist(lapply(test, function(x) x$cancertype%>% unique))
  test <- lapply(test, function(x) x%>%column_to_rownames("rowname")%>%dplyr::select(-cancertype)%>%as.data.frame)
  names(test) <- names
  gex_ccle <- test
  save(gex_ccle, file = "data/CCLE/CCLE_RNAseq_processed.RData")
}
#######################################


types <- names(table(cosmic$tissue)[table(cosmic$tissue) > 15 & table(cosmic$tissue) <100])
which.cancer.type <- script_params(vector = types, submission = T, args = args, parallel = F)

print(which.cancer.type)
msamples <- read.csv(paste0(path,"SraRunTable.txt"))
ccle <- init_annotations_preprocess(path.to.drugs_original = path,
                                    path.to.cosmic = "data/cosmic.RData")
msamples$name <- toupper(gsub(" ","",gsub("[[:punct:]]", "", msamples$Cell_Line)))
ccle$cl_info$meth <- lapply(ccle$cl_info$ccl_name, function(x) paste0("",as.character(msamples$Run[msamples$name == x & !is.na(msamples$name)]))) %>% unlist

# get samples only not in gdsc
not_in_gdsc <- row.names(ccle$resp[[which.cancer.type]])[lapply(row.names(ccle$resp[[which.cancer.type]]), function(x) ccle$cl_info$tissue_gdsc[ccle$cl_info$master_ccl_id == x] != which.cancer.type)%>%as.character == "TRUE"]
not_in_gdsc <- unlist(lapply(not_in_gdsc, function(x) ccle$cl_info$ccl_name[ccle$cl_info$master_ccl_id == x]%>%as.character))
samples_not_in_gdsc <- msamples$Run[which(msamples$name %in% not_in_gdsc)]%>%as.character; #cat(samples) # only samples not in gdsc


### Select all methylation samples for all cancer types cancertypes
samples <- lapply(
  types,
  function(type) msamples$Run[which(msamples$name %in% unlist(lapply(row.names(ccle$resp[[type]]), function(x) ccle$cl_info$ccl_name[ccle$cl_info$master_ccl_id == x]%>%as.character)))]%>%as.character
); names(samples)=types
#cat(samples$NB)
#cat(unlist(samples[!names(samples) %in% c("OV","NB")])[1:100])

### Import the methylation data
cancer.type <- which.cancer.type
library(methylKit, lib.loc = "~/R/x86_64-redhat-linux-gnu-library/3.6")
sample_list <- list.files("data/RRBS/", pattern = "sorted.bam", full.names = T)
sample_list <- sample_list[unlist(lapply(1:length(sample_list),function(i) any(unlist(lapply(samples[[cancer.type]], function(x) grepl(x,sample_list[i],fixed=T))))))] # only for special cancer type
sample_list <- sample_list[file.info(sample_list)$size > 100] # filter for size
print(sample_list)
dir.create(path = paste0("data/methylDB/",which.cancer.type))
path.to.save <- paste0("data/methylDB/",cancer.type,"/",cancer.type,"_ccle_meth.RData") # <- for counting ccle cell lines
if(file.exists(path.to.save)){
  load(path.to.save)
}else{
  my.methRaw=processBismarkAln(location = lapply(sample_list, function(x)x),
                               sample.id=lapply(lapply(unlist(lapply(sample_list, function(y) strsplit(y,"_")[[1]][1])),function(z) strsplit(z,"//")[[1]][2]), function(x)x),
                               assembly="hg19", 
                               read.context="CpG",
                               save.folder = paste0("data/methylDB/",which.cancer.type,"/"),
                               mincov = 3,
                               treatment = rep(1, length(sample_list))
  )
  save(my.methRaw, file=path.to.save)
}

#getMethylationStats(my.methRaw[[4]],plot=TRUE,both.strands=TRUE)
#getCoverageStats(my.methRaw[[2]],plot=TRUE,both.strands=TRUE)

meth_ccle <- methylKit::unite(my.methRaw, min.per.group = 6L)
meth_ccle@.Data[[1]] <- unlist(lapply(meth_ccle@.Data[[1]], function(x) chr_names$V2[chr_names$V1 == x])) # lift chromosomes
tmp <- meth_ccle@.Data[[1]]
if(any(duplicated(paste(tmp, meth_ccle@.Data[[2]],sep=".")))){
  tmp[duplicated(paste(tmp, meth_ccle@.Data[[2]],sep="."))] <- as.character(1:as.numeric(table(duplicated(paste(tmp, meth_ccle@.Data[[2]],sep=".")))["TRUE"])) # manually
}
meth_ccle@.Data[[1]] <- tmp

meth_ccle_values  <- percMethylation(meth_ccle, rowids = T)/100 # logit-transform beta-values
beta.val <- meth_ccle_values
Problems <- which(beta.val < 0 | beta.val > 1)
beta.val <- as.matrix(beta.val)
if (length(Problems) != 0) {
  beta.val[Problems] <- NA
}
onevalues <- which(beta.val == 1)
zerovalues <- which(beta.val == 0)
if (length(onevalues) > 0 | length(zerovalues) > 0) {
  if (length(onevalues) > 0) {
    beta.val[onevalues] <- NA
    beta.val[onevalues] <- max(beta.val, na.rm = T)
  }
  if (length(zerovalues) > 0) {
    beta.val[zerovalues] <- NA
    beta.val[zerovalues] <- min(beta.val, na.rm = T)
  }
}
beta.val <- log(beta.val/(1 - beta.val))
meth_ccle_values <- beta.val; rm(beta.val)
meth_ccle_granges <- as(meth_ccle,"GRanges")

df$row <- 1:nrow(df)
df_granges <- makeGRangesFromDataFrame(df, keep.extra.columns = T)
overlaps <- findOverlaps(query = df_granges,subject = meth_ccle_granges) # no filter for cancer type


### Analysis per cancer type
validation_cancertype <- cancer.type
S <- get_data_for_scatter(df = df, results_collected = results_collected, which.cancer.type = validation_cancertype)
df_overlap <- makeGRangesFromDataFrame((df[unique(overlaps@from),]) %>% dplyr::filter(cancertype == validation_cancertype), keep.extra.columns = T)
overlap <- findOverlaps(query = df_overlap,subject = meth_ccle_granges) # only for overlaps

#Import gex 
load("data/CCLE/CCLE_RNAseq_processed.RData")
ccle_gex <- gex_ccle[[cancer.type]]
row.names(ccle_gex) <- unlist(lapply(row.names(ccle_gex), function(x) paste0("",as.character(ccle$cl_info$master_ccl_id[ccle$cl_info$meth == x])))) #use ccle manifest for mapping cl IDs

#Import TCGA data
if(F){
  tcga_gex <- loadRData(paste(METADATA_TCGA, which.cancer.type,"_gex.RData",sep=""))
  candidate.rownames <- unlist(lapply(row.names(tcga_gex), function(x) paste0(strsplit(x,"-")[[1]][1:4], collapse = "-")))
  tcga_gex <- tcga_gex[!duplicated(candidate.rownames),]
  row.names(tcga_gex) = candidate.rownames[!duplicated(candidate.rownames)]
  tcga_meth <- tryCatch(loadRData(paste(METADATA_TCGA, which.cancer.type,"_betas.RData",sep="")),error = function(e) NA)
  if(class(tcga_meth) != "logical"){
    tcga_meth <- tcga_meth$beta[tcga_meth$phenotype == "TP",]%>%t
    row.names(tcga_meth) <- unlist(lapply(row.names(tcga_meth), function(x) paste0(strsplit(x,"-")[[1]][1:4], collapse = "-")))
  }
  int <- intersect(row.names(tcga_meth),row.names(tcga_gex))
  tcga_meth <- tryCatch(tcga_meth[int,], error = function(e) NA)
  tcga_gex <- tcga_gex[int,]
  if(class(tcga_meth) == "logical"){
    tcga_meth <- tcga_gex
  }
}


#lapply(1:length(df_overlap),
#       function(i){
ccle_names <- tolower(gsub(" ","",str_replace_all(ccle$drugs_info$DRUG_NAME, "[[:punct:]]", "")))
gdsc_names <- tolower(gsub(" ","",str_replace_all(df_overlap$drug, "[[:punct:]]", "")))

ccle_data <- list();for(i in 1:length(df_overlap)){
  tmp_ccle_sites <- meth_ccle_granges[overlap@to[overlap@from == i],];
  tmp_ccle_sites <- paste(tmp_ccle_sites@seqnames, tmp_ccle_sites@ranges@start, tmp_ccle_sites@ranges@start,sep=".")
  tmp_meth <- t(meth_ccle_values[tmp_ccle_sites,,drop = F])
  row.names(tmp_meth) = unlist(lapply(row.names(tmp_meth), function(x) paste0("",as.character(ccle$cl_info$master_ccl_id[ccle$cl_info$meth == x])))) #use ccle manifest for mapping cl IDs

  if(gdsc_names[i] %in% ccle_names){ # if drug is in CCLE
    ccle_drug_id <- ccle$drugs_info$master_cpd_id[ccle_names == gdsc_names[i]]
    tmp_resp <- ccle$resp[[validation_cancertype]][,as.character(ccle_drug_id), drop = F]
    
    int_meth_gex <- intersect(row.names(tmp_meth), row.names(ccle_gex))
    tmp_gex <- ccle_gex[int_meth_gex,]
    tmp_meth_gex <- tmp_meth[int_meth_gex,,drop = F]
    
    int <- intersect(row.names(tmp_meth), row.names(tmp_resp))
    tmp_meth <- tmp_meth[int,,drop=F]
    tmp_resp <- tmp_resp[int,,drop=F]
    
    validation <- unlist(lapply(1:ncol(tmp_meth), function(x) if(nrow(na.omit(cbind(tmp_resp[,1, drop = F],tmp_meth[,x, drop = F])))>3){summary(lm((tmp_resp[,1]/16)~ tmp_meth[,x],))$coefficients[,"Estimate"][2]}else{NA})) %>% na.omit %>% mean # /16 cause of concentration range ###cor(tmp_meth[,x],tmp_resp/16,use = "complete.obs")###
    discovery <- df_overlap$effectsize[i]
    
    int_not_gdsc <- unlist(lapply(int, function(x) ccle$cl_info$meth[ccle$cl_info$master_ccl_id == x]))
    int_not_gdsc <- int_not_gdsc[int_not_gdsc %in% samples_not_in_gdsc]
    int_not_gdsc <- as.character(unlist(lapply(int_not_gdsc, function(x) ccle$cl_info$master_ccl_id[ccle$cl_info$meth == x])))
    tmp_meth_not_gdsc <- tmp_meth[int_not_gdsc,,drop=F]
    tmp_resp_not_gdsc  <- tmp_resp[int_not_gdsc,,drop=F]
    validation_not_gdsc <- unlist(lapply(1:ncol(tmp_meth_not_gdsc), function(x) if(nrow(na.omit(cbind(tmp_resp_not_gdsc[,1, drop = F],tmp_meth_not_gdsc[,x, drop = F])))>3){summary(lm((tmp_resp_not_gdsc[,1]/16)~ tmp_meth_not_gdsc[,x],))$coefficients[,"Estimate"][2]}else{NA})) %>% na.omit %>% mean
    
  
    
    tmp <- cbind(tmp_resp, tmp_meth)
    tmp$not_gdsc <- row.names(tmp) %in% int_not_gdsc
    tmp <- list(df = tmp, 
                discovery = discovery, 
                validation = validation,
                validation_not_gdsc = validation_not_gdsc,
                drug = as.character(df_overlap$drug[i]),
                index.in.df = df_overlap$row[i]#,
                #pattern$mg,
                #pattern_tcga$mg
                )
  }else{
    message(paste("Drug not found:",gdsc_names[i]))
    int_meth_gex <- intersect(row.names(tmp_meth), row.names(ccle_gex))
    tmp_gex <- ccle_gex[int_meth_gex,]
    tmp_meth_gex <- tmp_meth[int_meth_gex,,drop = F]

    
    tmp <-  list(df = NA, 
                 discovery = NA, 
                 validation = NA,
                 validation_not_gdsc = NA,
                 drug = NA,
                 index.in.df = df_overlap$row[i]#,
                 #pattern$mg,
                 #pattern_tcga$mg
    )
  }
  print(i)
  #return()
  ccle_data[[i]] <- tmp
}; names(ccle_data) = df_overlap$gene
#)
test <- lapply(ccle_data, function(x){
  tmp <- as.data.frame(x); if(class(x$df)!="logical"){names <- paste0("df.",colnames(x$df)); nest(tmp, names)}else{cbind(tmp[,-1],data = NA)}
  })
test <- do.call(bind_rows, test)

#ccle_data_summary <- do.call(rbind,lapply(ccle_data, function(obj) unlist(c(obj["discovery"], obj["validation"],obj["drug"]))))
ccle_data_summary <- test
save(ccle_data_summary, file=paste0("methylDB/",cancer.type,"/",cancer.type,"_ccle_meth_res_updated_3.RData"))






