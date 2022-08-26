## load packages required for analysis
library(limma)
library(minfi)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylation450kmanifest)
library(ChAMP)
library(readxl)

ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
RGset <- read.metharray.exp("data/GSE68379_RAW/", verbose=T) # read the test idat files from the test folder
xReactiveProbes <- read_excel("/metadata/48639-non-specific-probes-Illumina450k.xlsx") # filter cross-reactive probes

sentrix <- unlist(lapply(RGset@colData@rownames, function(x) paste(strsplit(x, "_")[[1]][-1],collapse="_")))
sample_sheet <- read_excel("data/methSampleId_2_cosmicIds.xlsx")
sample_sheet$meth <- unlist(lapply(1:nrow(sample_sheet), function(x) paste(sample_sheet$Sentrix_ID[x],sample_sheet$Sentrix_Position[x], sep="_")))

cosmics <- unlist(lapply(sentrix, function(x) sample_sheet$cosmic_id[x == sample_sheet$meth]))
stopifnot(length(cosmics) == 1028)
cancer_types <- lapply(cosmics, function(x) unlist(na.omit(unique(cosmic$tissue[cosmic$cosmic_id == x]))))
all_cancer_types <- names(table(cancer_types%>% unlist))


detP <- detectionP(RGset) # get detection pvalues
#barplot(colMeans(detP), las=2, cex.names=0.8, main="Mean detection p-values", ylim=c(0,0.051))  # look at the mean detection P-values across all samples to identify any failed samples
#abline(h=0.05,col="red")
keep <- colMeans(detP) < 0.05 # remove too high detection pvalue
RGset <- RGset[,keep]
detP <- detP[,keep]
GRset.funnorm <- preprocessNoob(RGset, verbose = T) # "noob normalization"
Mset.funnorm <- mapToGenome(GRset.funnorm) # map to genome build
#densityPlot(getBeta(Mset.funnorm), main="Funnorm") # beta density plot
#plotMDS(getM(Mset.funnorm), top=1000, gene.selection = 'common') # limma plot on principal components
detP <- detP[match(featureNames(Mset.funnorm),rownames(detP)),] # ensure that probes are in same order than Mset
keep <- rowSums(detP < 0.01) == ncol(Mset.funnorm)
Mset.funnorm.flt <- Mset.funnorm[keep,]
keep <- !(featureNames(Mset.funnorm.flt) %in% ann450k$Name[ann450k$chr %in% c("chrX","chrY")]) # filter out sex chromosomes
Mset.funnorm.flt <- Mset.funnorm.flt[keep,]
keep <- !(featureNames(Mset.funnorm.flt) %in% xReactiveProbes$TargetID)
Mset.funnorm.flt <- Mset.funnorm.flt[keep,]

myNorm <- champ.norm(mset = Mset.funnorm.flt, beta = getBeta(Mset.funnorm.flt), arraytype = "450K", cores = 1, method = "BMIQ")

betas <- myNorm%>%t
row.names(betas) = cosmics
stopifnot(nrow(betas) == 1028)

for(which in all_cancer_types){
  betas_type <- betas[as.character(cancer_types) == which,]
  save(betas_type , file = paste("metadata/","methyl_processed_",as.character(which),".RData",sep=""))
}
