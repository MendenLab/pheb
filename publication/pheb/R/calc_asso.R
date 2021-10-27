suppressPackageStartupMessages({
library(CpGassoc)
library(ggplot2)
library(ggrepel)
library(stringr)
})

mutate_DAT <- function(
### mutates the GDSC data object into aggregated beta-value matrix
######################################################
  dmr.object = NULL,
  data.object = NULL,
  annotation.object = NULL
){
  dmr.list <- dmr.object$DMR[[1]]
  regions <- lapply(1:nrow(dmr.list), function(x) unlist(dmr.list$start[x]:dmr.list$end[x]))
  regions_chrom <- dmr.list$seqnames
  cpgs <- lapply(1:length(regions), function(x) unlist(row.names(hm450.manifest.hg19)[(hm450.manifest.hg19$CpG_beg %in% regions[[x]] | hm450.manifest.hg19$CpG_end %in% regions[[x]]) & 
                                                                    hm450.manifest.hg19$CpG_chrm  ==  as.character(regions_chrom[x])        ]))
  print(fot(unlist(cpgs) %in% colnames(data.object$meth)))
  
  # Create the m-value dataframe
  object <- matrix(nrow = nrow(data.object$meth), ncol = nrow(dmr.list))
  colnames(object) = unlist(lapply(1:nrow(dmr.list), function(x) paste(dmr.list$seqnames[x],":",as.character(dmr.list$start[x]),"-",as.character(dmr.list$end[x]), sep = "")))
  row.names(object) = row.names(data.object$meth)
    
  for(i in 1:ncol(object)){
    temp <- data.object$meth[, colnames(data.object$meth) %in% cpgs[[i]]]
    temp <- apply(temp, 1, mean)
    object[,i] <- temp
  }
  
  gene_anno <- lapply(cpgs, function(x) annotation.object$reference_features$UCSC_RefGene_Name[annotation.object$reference_features$IlmnID %in% unlist(x)] )
  data.object$method <- dmr.object$method
  data.object$meth <- object
  data.object$gene_anno <- gene_anno
  
  return(data.object)
}######################################################


asso <- function(
  ### returns dataframe of single association test
  ######################################################
  data.object = NULL,
  annotation.object = NULL,
  which.drugs.index = T, #
  metadata = F, # must be =F
  description = "", # descriptional string for the uniqueness of the metadata object
  min_screened = 4, # minimal cell lines screened per drug for testing
  min_auc = 0.7, # auc threshold for being a responder
  min_responder = 3, # how many responder minimally
  logit.transform = T,
  scaling = F,
  covariates = NULL,
  meth = T, # only =T, because Gex analysis is deprecated
  lrtest = F #if not false, need to supply names in covariates for lrtest
){
  if(!all(unlist(lapply(which.drugs.index, function(x) is.integer(x)))))
    which.drugs.index <- 1:ncol(data.object$resp)
  which.drugs.index <- na.omit(which.drugs.index[(apply(data.object$resp,2 , function(x) 
    length(which(!is.na(x))) > min_screened &
    length(which(x < min_auc)) >= min_responder )) 
    ]) # get the drugs with cell lines screened more than min_screened
  
  if(!(lrtest == F)){
    cpg.assoc_my <- dget("R/cpg.assoc_my_lrtest.R") # assign the function such that it executes LR-tests instead of single variate
    print("Exectute LR-Test !")
  }
  
  if(meth){
    resp <- data.object$resp[intersect(row.names(data.object$resp),row.names(data.object$meth)),,drop=F]
    if(all(unlist(lapply(which.drugs.index, function(x) is.integer(x))))) 
      resp <- resp[,which.drugs.index, drop = F]
    obs <- t(data.object$meth[intersect(row.names(data.object$resp),row.names(data.object$meth)),])
    results <- NULL
    gene_names <- lapply(data.object$gene_anno, function(x) unique(as.character(x)))
    gene_names <- lapply(gene_names, function(x) str_c(x, sep = ";", collapse = ";"))
    gene_names <- lapply(gene_names, function(x) as.character(str_c(unique(unlist(strsplit(x,";"))),collapse = " ")))
    
    print(paste("Performing on number of drugs: ",ncol(resp), sep=""))
    if(length(resp)!=0){ #!
      if(ncol(resp)!=0){ #!
        obs <- t(obs) # was in loop before
        if(!is.null(covariates)){covariates <- covariates[row.names(obs),,drop=F]}else{print("No Covariates supplied !")}
        for(i in 1:ncol(resp)){
          if(metadata == FALSE | !file.exists(paste("metadata/cpg_asso/",data.object$which,"_",as.character(colnames(resp)[i]),"_",description,".RData",sep=""))){
            if(scaling == T){
              PHENOTYPE <- scale(resp[,i])
            } else {
              PHENOTYPE <- resp[,i]
            }
            
            names(PHENOTYPE) <- row.names(resp) ## new
            temp <- cpg.assoc_my(beta.val = obs, pheno = PHENOTYPE, logit.transform = logit.transform, covariates = covariates, lrtest = lrtest)
            #save(temp, file = paste("metadata/cpg_asso/",data.object$which,"_",data.object$DMR$method,"_",as.character(colnames(resp)[i]),"_",description,".RData",sep=""))
          } else {
            load(paste("metadata/cpg_asso/",data.object$which,"_",data.object$DMR$method,"_",as.character(colnames(resp)[i]),"_",description,".RData",sep=""))
          }
          
          effectsize <- temp$effectsize
          temp$effectsize <- effectsize
          temp$indep <- rep(colnames(resp)[i], nrow(temp)) # annotate my drug ID
          drug_name <- as.character(annotation.object$drugs_info$DRUG_NAME[annotation.object$drugs_info$DRUG_ID == strsplit(colnames(data.object$resp)[which.drugs.index[i]],"-")[[1]][1]])
          temp$drug <- rep(drug_name, nrow(temp)) # annotate my drug name
          temp$gene <- if(is.list(data.object$gene_anno)){gene_names}else{temp$CPG.Labels} # annotate genes if it is suppliedwarnings
          amount <- length(which(!is.na(resp[,i]))) # how many drugs are not NA
          temp$amount <- rep(amount, nrow(temp))
          temp$confidence <- -log2(unlist(apply(obs, 2, function(x) sd(x, na.rm = T))))# - log2(sd_BT(resp[,i])) #confidence measure, can be anything
          # changed 1 to 2, since obs is now transposed
          results <- rbind(results, temp)
          print(paste("Progress: ", as.character(round(i/ncol(resp)*100)),"%", sep=""))
        }
      }
    } else {results <- NULL}
  }
  return(results)
}######################################################


asso_mut_meth <- function(
  ### returns dataframe of single association test
  ######################################################
  data.object = NULL,
  annotation.object = NULL,
  which.drugs.index = T, # 
  metadata = F, # must be =F
  description = "", # descriptional string for the uniqueness of the metadata object
  min_screened = 4, # minimal cell lines screened per drug for testing
  min_auc = 0.7, # auc threshold for being a responder
  min_responder = 3, # how many responder minimally
  logit.transform = T,
  scaling = F,
  covariates = NULL,
  meth = T, # only =T, because Gex analysis is deprecated
  lrtest = F # if not false, need to supply names in covariates for lrtest
){
  if(!all(unlist(lapply(which.drugs.index, function(x) is.integer(x))))) 
    which.drugs.index <- 1:ncol(data.object$resp)
  which.drugs.index <- na.omit(which.drugs.index[(apply(data.object$resp,2 , function(x) 
    length(which(!is.na(x))) > min_screened &
      length(which(x < min_auc)) >= min_responder )) 
    ]) # get the drugs with cell lines screened more than min_screened
  
  if(!(lrtest == F)){
    cpg.assoc_my <- dget("R/cpg.assoc_my_lrtest.R") # assign the function such that it executes LR-tests instead of single variate
    print("Exectute LR-Test !")
  }
  
  if(meth){
    resp <- data.object$resp[intersect(row.names(data.object$resp),row.names(data.object$meth)),,drop=F]
    if(all(unlist(lapply(which.drugs.index, function(x) is.integer(x))))) # 
      resp <- resp[,which.drugs.index, drop = F]
    obs <- t(data.object$meth[intersect(row.names(data.object$resp),row.names(data.object$meth)),])
    results <- NULL
    gene_names <- lapply(data.object$gene_anno, function(x) unique(as.character(x)))
    gene_names <- lapply(gene_names, function(x) str_c(x, sep = ";", collapse = ";"))
    gene_names <- lapply(gene_names, function(x) as.character(str_c(unique(unlist(strsplit(x,";"))),collapse = " ")))
    
    print(paste("Performing on number of alterations: ",ncol(resp), sep=""))
    if(length(resp)!=0){ #!
      if(ncol(resp)!=0){ #!
        obs <- t(obs) # was in loop before
        if(!is.null(covariates)){covariates <- covariates[row.names(obs),,drop=F]}else{print("No Covariates supplied !")}
        for(i in 1:ncol(resp)){
          if(metadata == FALSE | !file.exists(paste("metadata/cpg_asso/",data.object$which,"_",as.character(colnames(resp)[i]),"_",description,".RData",sep=""))){
            if(scaling == T){
              PHENOTYPE <- scale(resp[,i])
            } else {
              PHENOTYPE <- resp[,i]
            }
            
            names(PHENOTYPE) <- row.names(resp) ## new
            temp <- assoc_mut_meth(beta.val = obs, pheno = PHENOTYPE, logit.transform = logit.transform, covariates = covariates, lrtest = lrtest)
            #save(temp, file = paste("metadata/cpg_asso/",data.object$which,"_",data.object$DMR$method,"_",as.character(colnames(resp)[i]),"_",description,".RData",sep=""))
          } else {
            load(paste("metadata/cpg_asso/",data.object$which,"_",data.object$DMR$method,"_",as.character(colnames(resp)[i]),"_",description,".RData",sep=""))
          }
          
          
          effectsize <- temp$effectsize
          temp$effectsize <- effectsize
          temp$indep <- rep(colnames(resp)[i], nrow(temp)) # annotate my drug ID
          #drug_name <- as.character(annotation.object$drugs_info$DRUG_NAME[annotation.object$drugs_info$DRUG_ID == strsplit(colnames(data.object$resp)[which.drugs.index[i]],"-")[[1]][1]])
          #temp$drug <- rep(drug_name, nrow(temp)) # annotate my drug name
          temp$gene <- if(is.list(data.object$gene_anno)){gene_names}else{temp$CPG.Labels} # annotate genes if it is suppliedwarnings
          amount <- length(which(!is.na(resp[,i]))) # how many drugs are not NA
          temp$amount <- rep(amount, nrow(temp))
          temp$confidence <- -log2(unlist(apply(obs, 2, function(x) sd(x, na.rm = T))))# - log2(sd_BT(resp[,i])) #confidence measure, can be anything
          # changed 1 to 2, since obs is now transposed
          results <- rbind(results, temp)
          #cat(paste("-", as.character(round(i/ncol(resp)*100)),"", sep=""))
        }
      }
    } else {results <- NULL}
  }
  return(results)
}######################################################


### PLOT FUNCTIONS

plot_volcano <- function(
### plots volcano plot
######################################################
  df, 
  Cond, 
  sign, 
  annotations='', 
  annotation.object = NULL, 
  pval = 0.05, 
  eff = 5
){
  
  drugs_info <- annotation.object$drugs_info
  df <- df[Cond,]
  significant <- sign[Cond]
  
  if(annotations=='')
    annotations = rep(F, length(sign))
  gene_annotation <- (df$P.Value <= sort(df$P.Value[significant], decreasing = F)[min(30,length(which(significant)))]) #for greatest significant hits
  df$gene_annotation <- gene_annotation
  df$significant <- significant
  p1 <- ggplot(df, aes(x=effectsize, 
                       y=-log10(P.Value),
                       label= ifelse( gene_annotation & significant,gene,'')
                       )) + 
    geom_point(aes(size=df$confidence, 
                   color = factor(annotation), 
                   shape = factor(annotation)
                   ), alpha = 0.5) + #size=df$amount; minus confidence measure (is -log(sd))
    scale_size(range = c(1, 8)) +
    #labs(title = paste("Significant drug-island association for ",annotation.object$which,sep=""), colour="Drug", size="Confidence score") +
    xlab("Effect size") +
    scale_shape_manual(values = c(((0:21)[-c(19,20,21)]),24,25,23))+
    #ylab("-log10(p)") + ylim(c(-0.,10)) + xlim(c(-125,125)) + 
    theme(plot.title = element_text(hjust = 0.5))+
    #if(length(which(gene_annotation & significant)) !=0){
      geom_text_repel(
        data = subset(df, effectsize <0),
        size = 2,
        hjust = 1,
        direction = "y",
        nudge_x = -1.5 - subset(df, effectsize<0)$effectsize,
        segment.size = 0.2
      )+geom_text_repel(
          data = subset(df, effectsize >0),
          size = 2,
          hjust = 0,
          direction = "y",
          nudge_x = 1.5 - subset(df, effectsize>0)$effectsize,
          segment.size = 0.2
        )+
    theme_minimal()+
    labs(shape = "Cancer type",col = "Cancer type", size = "# of probes")
  return(p1)
}######################################################


plot_volcano_2 <- function(
### plots volcano plot
######################################################
df, Cond, annotations='', annotation.object = NULL, pval = 0.05, eff = 5,
max_p = 8.5, max_eff = 120, name = ""
){
  drugs_info <- annotation.object$drugs_info
  df <- df[Cond,]
  p2 <- p1 + geom_point(data=df, aes(x=effectsize, y=-log10(P.Value), size=df$confidence, label=""), shape=16, alpha = 0.05) + # negaive confidence measure
    #geom_vline(xintercept = eff, linetype="dashed",color = "black", size=0.3)+
    #geom_vline(xintercept = -eff, linetype="dashed",color = "black", size=0.3)+
    #geom_hline(yintercept = -log10(pval), linetype="dashed",color = "black", size=0.3)+
    labs(title = paste(annotation.object$which,": Significant drug-island associations for ","",sep=""), colour="Drug", size="Confidence score") +
    xlab("Effect Size") +
    ylab("-log10(p)") + ylim(c(-0.,max_p)) + xlim(c(-max_eff,max_eff)) + 
    theme(plot.title = element_text(hjust = 0.5))
  return(p2)
}


### Basically just switched the indep with dep
assoc_mut_meth <- function(beta.val = RESULTS$DAT_tcga$meth, pheno = RESULTS$DAT_tcga$resp[,400], covariates = NULL, logit.transform = F, logit.transform.pheno = F, 
                            lrtest = F) 
{
  betaval <- as.matrix(beta.val)
  overlapping <- intersect(row.names(beta.val), names(pheno))
  beta.val <- betaval[overlapping, ]
  pheno <- pheno[overlapping]
  if (!is.null(covariates)) {
    covariates <- covariates[overlapping, , drop = F]
  }
  if (logit.transform) {
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
  }
  if (logit.transform.pheno) {
    Problems <- which(pheno < 0 | pheno > 1)
    if (length(Problems) != 0) {
      pheno[Problems] <- NA
    }
    onevalues <- which(pheno == 1)
    zerovalues <- which(pheno == 0)
    if (length(onevalues) > 0 | length(zerovalues) > 0) {
      if (length(onevalues) > 0) {
        pheno[onevalues] <- NA
        pheno[onevalues] <- max(pheno, na.rm = T)
      }
      if (length(zerovalues) > 0) {
        pheno[zerovalues] <- NA
        pheno[zerovalues] <- min(pheno, na.rm = T)
      }
    }
    pheno <- log(pheno/(1 - pheno))
  }
  covariates.not.na <- covariates[!is.na(pheno), , drop = F]
  covariates.not.na <- apply(covariates.not.na, 2, function(x) length(levels(factor(x))) > 
                               1)
  covariates <- covariates[, covariates.not.na, drop = F]
  df <- data.frame(matrix(ncol = 3, nrow = 0))
  colnames(df) <- c("P.Value", "effectsize", "t-value")
  for (i in 1:ncol(beta.val)) {
    if (is.null(covariates)) {
      lm2 <- tryCatch(lm(betaval[, i] ~ factor(pheno)), error = function(e) NULL)
    }
    else {
      lm2 <- tryCatch(lm(betaval[, i] ~ factor(pheno) + ., data = covariates), 
                      error = function(e) NULL)
    }
    lm2 <- summary(lm2)
    #print(lm2)
    if (length(which(!is.na(betaval[, i]))) > 2 & min(table(factor(pheno)))>3 & length(table(factor(pheno)))==2) {
      df[i, ] <- c(lm2$coefficients[2, 4], lm2$coefficients[2, 
                                                            1], lm2$coefficients[2, 3])
    }
    else {
      df[i, ] <- c(NA, NA, NA)
    }
  }
  fdr <- p.adjust(df$P.Value, "BH")
  df$FDR <- fdr
  df$CPG.Labels <- colnames(betaval)
  return(df)
}




