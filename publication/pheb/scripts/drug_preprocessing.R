source("R/script_params.R")
load(paste("data/cosmic.RData","",sep=''))
types <- names(table(cosmic$tissue)[table(cosmic$tissue) > 15 & table(cosmic$tissue) <100])
which.cancer.type <- script_params(vector = "BRCA", # any cancer type works
                                   submission = T, args = 1, parallel = F)

# Create annotation objects for GDSC drugs and methylation array
ANN <- init_annotations(which.cancer.type,
                        drug=F,
                        meth=F,
                        cfe=F, 
                        path.to.cosmic= COSMIC_PATH, 
                        path.to.drugs_original2 = GDSC2_SCREEN_PATH,
                        path.to.drugs_original1= GDSC1_SCREEN_PATH,
                        path.to.drugs_info = DRUG_INFO_PATH,
                        path.to.reference_original = METH_MANIFEST_PATH)



cancer_types <- as.data.frame(sort(table(ANN$cosmic['tissue']), decreasing = TRUE))
names_type <- function(can_ty){
  nam <- ANN$cosmic[ANN$cosmic$tissue==can_ty & !is.na(ANN$cosmic$tissue),"cosmic_id"]
  return(nam)}

### Import drug response ####
#drugs_original <- as.data.frame(read.csv("data/GDSC_drugs/PANCANCER_IC_Thu Aug  1 12_17_03 2019.csv"))
#drugs_original_original <- as.data.frame(read.csv("data/GDSC_drugs/PANCANCER_IC_Thu Aug  1 13_41_49 2019.csv"))
drugs_original <- read_excel("data/GDSC_drugs/GDSC2_fitted_dose_response_17Jul19.xlsx")
drugs_original_original <- read_excel("data/GDSC_drugs/GDSC1_fitted_dose_response_17Jul19.xlsx")
########################



### Preprocess ####
drugs_2 <- drugs_original[,c("COSMIC_ID", "CELL_LINE_NAME", "DRUG_ID","DRUG_NAME","MAX_CONC","LN_IC50" ,"AUC", "DATASET")]
drugs_2_original <- drugs_original_original[,c("COSMIC_ID", "CELL_LINE_NAME", "DRUG_ID","DRUG_NAME","MAX_CONC","LN_IC50" ,"AUC", "DATASET")]
drugs_2 <- rbind(drugs_2_original, drugs_2)
colnames(drugs_2) = c("COSMIC_ID","CELL_LINE_NAME","DRUG_ID","DRUG_NAME","MAX_CONC_MICROMOLAR","LN_IC50","AUC","Dataset.version")
drugs_2$CELL_LINE_NAME <- make.names(drugs_2$CELL_LINE_NAME)
drugs_2$DRUG_ID <- as.factor(drugs_2$DRUG_ID)
drugs_2$COSMIC_ID <- as.character(drugs_2$COSMIC_ID)
drugs_2$"MY_ID" <- factor(paste(as.character(drugs_2$DRUG_ID), as.character(drugs_2$Dataset.version), sep="-"))
drugs_2 <- drugs_2[c("COSMIC_ID","MY_ID","AUC")]
drugs_spread <- spread(drugs_2, COSMIC_ID, AUC)
#row.names(drugs_spread) = drugs_spread$COSMIC_ID
#drugs_spread <- drugs_spread[-1]
###################


### Create obverables for cell lines ####
FEATURE <- "MY_ID"
for(i in 1:nrow(cancer_types)){
  response <- as.data.frame(drugs_spread[FEATURE])
  colnames(response) = FEATURE
  for(j in 1:length(names_type(cancer_types[i,1]))){
    
    if(names_type(cancer_types[i,1])[j] %in% colnames(drugs_spread)){
      tmp <- as.data.frame(drugs_spread[,names_type(cancer_types[i,1])[j]])
      colnames(tmp) = names_type(cancer_types[i,1])[j]
      response <- cbind(response, tmp)
      print(names_type(cancer_types[i,1])[j])
    }
  }
  row.names(response) = response[,FEATURE]
  response <- response[-1] # moved the FEATURE column to the rownames
  response <- as.data.frame(t(response)) # transpose such that cell lines are rows
  filename <- paste(paste("metadata/GDSC_create_drugs_AUC_2/",paste("response",FEATURE,sep = '_'), sep=''),"_",cancer_types[i,1],".RData", sep='')
  print(filename)
  save(response, file = filename)
  if(i==1){pancancer <- response} else {pancancer <- rbind(pancancer, response)}
}
filename <- paste(paste("metadata/GDSC_create_drugs_AUC_2/",paste("response",FEATURE,sep = '_'), sep=''),"_","PANCANCER",".RData", sep='')
print(filename)
save(pancancer, file = filename)
##############################################