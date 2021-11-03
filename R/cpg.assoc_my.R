cpg.assoc_my <- function(
  beta.val = RESULTS$DAT_tcga$meth, #data.frame
  pheno = RESULTS$DAT_tcga$resp[,400], #named vector
  covariates = NULL, # RESULTS$DAT$covariates
  logit.transform = F,
  logit.transform.pheno = F,
  lrtest = F
){
  betaval <- as.matrix(beta.val)
  overlapping <- intersect(row.names(beta.val),names(pheno))
  beta.val <- betaval[overlapping,]
  pheno <- pheno[overlapping]
  if(!is.null(covariates)){covariates <- covariates[overlapping,, drop= F]}
  
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
    #beta.val <- as.matrix(pheno)
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
  
  covariates.not.na <- covariates[!is.na(pheno),,drop = F]
  covariates.not.na <- apply(covariates.not.na,2, function(x) length(levels(factor(x)))>1)
  covariates <- covariates[,covariates.not.na,drop=F]
  if(ncol(covariates)== 0){
    covariates <- NULL # added to protect against full removal
  }
  
  df <- data.frame(matrix(ncol = 3, nrow = 0))
  colnames(df) <- c("P.Value","effectsize","t-value")
  for(i in 1:ncol(beta.val)){ # the main loop with takes up lots of computational time
    #print(i)
    if(is.null(covariates)){
      lm2 <- tryCatch(lm(pheno ~ betaval[,i]), error = function(e) NULL)
    }else{
      lm2 <- tryCatch(lm(pheno ~ betaval[,i] + ., data = covariates), error = function(e) NULL)
    }
    lm2 <- summary(lm2)
    #print(lm2)
    if(length(which(!is.na(betaval[,i])))>2){ 
      df[i,] <- c(lm2$coefficients[2,4], lm2$coefficients[2,1], lm2$coefficients[2,3])
      #print(c(lm2$coefficients[2,4], lm2$coefficients[2,1], lm2$coefficients[2,3]))
    } else {
      df[i,] <- c(NA,NA,NA)
    }
  }
  fdr <- p.adjust(df$P.Value, "BH")
  df$FDR <- fdr
  df$CPG.Labels <- colnames(betaval)
  
  return(df)
}
