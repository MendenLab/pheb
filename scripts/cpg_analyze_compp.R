library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
data(Locations)

prepare_csv <- function(
  data.Locations = Locations,
  files.list = files.list,
  replace = T,
  which.cancer.type = "SKCM"
){
  print(paste0("Number of files to mine: ",as.character(length(files.list))))
  for(files.list.index in 1:length(files.list)){
    test <- loadRData(files.list[files.list.index])
    if(files.list.index == 1){
      tmp <- test$DAT$resp
    }else{ # make the response dataframe <------
      #tmp <- cbind(tmp,test$DAT$resp)
    }
    drugs.list <- unique(test$cpg_asso$indep)
    print(paste0("Number of drugs: ",as.character(length(drugs.list))))
  
    if(replace == F){
      if(length(drugs.list)>0){
        logic <- as.logical(lapply(1:length(drugs.list), 
                    function(x) file.exists(paste0(dirname(files.list[files.list.index]),"/",which.cancer.type,".meth.",as.character(drugs.list[x]),".csv")) ))
        print(paste0(as.character(length(which(logic))),"/",as.character(length(logic))," exists !"))
      }else{
        print(paste0("0","/","0"," exists !"))
        logic <- FALSE
      }
      if(length(logic)!= length(drugs.list) & length(drugs.list)!=0)
        stop("Something is wrong !")
      
      drugs.list <- drugs.list[!logic]
    }
    if(length(drugs.list >0 )){
      for(drugs.list.index in 1:length(drugs.list)){
        cpgs <- test$cpg_asso$CPG.Labels[test$cpg_asso$indep == drugs.list[drugs.list.index]]
      
        Locations.local <- Locations[cpgs,]
        Locations.local <- data.frame(chrom=Locations.local$chr, 
                              start=Locations.local$pos, 
                              end=Locations.local$pos,
                              name = test$cpg_asso$CPG.Labels[test$cpg_asso$indep == drugs.list[drugs.list.index]],
                              pvalue=test$cpg_asso$P.Value[test$cpg_asso$indep == drugs.list[drugs.list.index]],
                              effectsize=test$cpg_asso$effectsize[test$cpg_asso$indep == drugs.list[drugs.list.index]]
                              )
        #Locations.local$pvalue <- p.adjust(Locations.local$pvalue, method = "BH")  #<- likely, combp adjusts for fdr by itself
        #colnames(Locations.local) = c("name","chrom", "start","end","pvalue")
        Locations.local <- Locations.local[,c("chrom", "start","end","name","pvalue","effectsize")]
        #ch <- unlist(lapply(Locations.local$chrom, function(x) gsub("chr","",x) ) )
        Locations.local <- Locations.local[order(Locations.local$chrom, Locations.local$start),]
      
        write.table(x = Locations.local, file = paste0(dirname(files.list[files.list.index]),"/",which.cancer.type,".meth.",as.character(drugs.list[drugs.list.index]),".csv"),
                    row.names = F, sep = "\t", dec = ".")
      }
    }
    progress <- round(files.list.index/length(files.list)*100)
    print(paste0("Progress: ",as.character(progress)," %"))
  }
  return(tmp)
}
response <- prepare_csv(data.Locations = Locations, 
                        files.list = files.list, 
                        replace = F, 
                        which.cancer.type = which.cancer.type)


### Perform the comb-p DMR search
system(paste("miniconda3/bin/conda init;
              . miniconda3/etc/profile.d/conda.sh;
              conda activate;
              cd R/;
             ./compp.sh",as.character(which.cancer.type)))


drugs.list <- list.files(path.to.metadata, pattern=paste0(which.cancer.type,"\\..+.csv.combp.regions.bed.gz"), all.files=FALSE,
                         full.names=TRUE)
drugs.list.names <- list.files(path.to.metadata, pattern=paste0(which.cancer.type,"\\..+.csv.combp.regions.bed.gz"), all.files=FALSE)
xtest <- length(drugs.list.names)
drugs.csv <- lapply(drugs.list,
                    function(x) gsub(".combp.regions.bed.gz", "", x)) %>% unlist
drugs.list.names <- lapply(drugs.list.names,
  function(x) gsub(".csv.combp.regions.bed.gz", "", gsub(paste0(as.character(which.cancer.type),"\\.meth\\."),"", x))) %>% unlist
if(length(drugs.list.names) != xtest)
  stop("Something is wrong !")
dmrs <- as.data.frame(matrix(ncol = 5, nrow=0))
for(i in 1:length(drugs.list)){
  temp <- NULL
  temp <- tryCatch(read.table(gzfile(drugs.list[i]), sep = "\t"), error= function(e) NULL)
  if(!is.null(temp)){
    temp$drug <- drugs.list.names[i]
  }
  if(!is.null(temp)){
    dmrs <- rbind(dmrs, temp)
  }
}

# Add annotations
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann450k = getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
colnames(dmrs) = c("chrom", "start","end","pvalue","nprobes","drug")
dmrs$index <- 1:nrow(dmrs)
dmrs$cpg <- lapply(1:nrow(dmrs), function(x)
  row.names(ann450k)[ann450k$pos %in% dmrs[x,]$start:dmrs[x,]$end & ann450k$chr == dmrs[x,]$chrom])
dmrs$enchancer <- lapply(1:nrow(dmrs), function(x)
  ann450k[unlist(dmrs$cpg[[x]]),"Enhancer"])
dmrs$dhs <- lapply(1:nrow(dmrs), function(x)
  ann450k[unlist(dmrs$cpg[[x]]),"DHS"])
dmrs$island <- lapply(1:nrow(dmrs), function(x)
  unique(ann450k[unlist(dmrs$cpg[[x]]),"Islands_Name"])
  )
dmrs$gene <- lapply(1:nrow(dmrs), function(x) 
  as.character(unlist(lapply(
  ann450k[unlist(dmrs$cpg[[x]]),"UCSC_RefGene_Name"],
    function(x) unique(unlist(strsplit(x, ";")[[1]]))
  )))
  )
dmrs$regulator <- lapply(1:nrow(dmrs), function(x)
  ann450k[unlist(dmrs$cpg[[x]]),"Regulatory_Feature_Group"])

dmrs$genepos <- lapply(1:nrow(dmrs), function(x)
  ann450k[unlist(dmrs$cpg[[x]]),"UCSC_RefGene_Group"])

probearray1 <- list()
probearray2 <- list()
for(j in unique(dmrs$drug)){
  assoo <- read.csv(drugs.csv[which(drugs.list.names %in% j)], sep = "\t")
  row.names(assoo) = assoo$name
  dmrscpg <- dmrs$cpg[dmrs$drug == j]
  pvalues <- lapply(dmrscpg, function(x) assoo[x,"pvalue"])
  effectsizes <- lapply(dmrscpg, function(x) assoo[x,"effectsize"])
  probearray1 <- c(probearray1, pvalues)
  probearray2 <- c(probearray2, effectsizes)
}

dmrs$pvalues <- probearray1
dmrs$effectsizes <- probearray2
dmrs$eff <- lapply(dmrs$effectsizes, function(x) as.numeric(na.omit(unique(sign(x)))))


# SAVE RESULTS
path.to.metadata <- "metadata/results_dmp_new"
save(dmrs, file = paste0(path.to.metadata,"_combp/",which.cancer.type,".regions.RData"))
save(response, file = paste0(path.to.metadata,"_combp/",which.cancer.type,".response.RData"))


