savepath <- "metadata"

source("R/script_params.R")
load(paste("metadata/cosmic.RData","",sep=''))
types <- names(table(cosmic$tissue)[table(cosmic$tissue) > 15 & table(cosmic$tissue) <100])


if(!file.exists(paste0(savepath, "/ctrp_validation.RData"))){
  which.cancer.type <- script_params(vector = types,
                                     submission = T, args = NA, parallel = F) #index for cancer type ! was args[2], number of containers
  ccle <- init_annotations_preprocess(path.to.drugs_original = path,
                                      path.to.cosmic = "data/cosmic.RData")
  msamples$name <- toupper(gsub(" ","",gsub("[[:punct:]]", "", msamples$Cell_Line)))
  ccle$cl_info$meth <- lapply(ccle$cl_info$ccl_name, function(x) paste0("",as.character(msamples$Run[msamples$name == x & !is.na(msamples$name)]))) %>% unlist
  
  types <- names(table(cosmic$tissue)[table(cosmic$tissue) > 15 & table(cosmic$tissue) <100])
  ctrptmp <- list(); for(type in types){
    print(type)
    ccle_names <- tolower(gsub(" ","",str_replace_all(ccle$drugs_info$DRUG_NAME, "[[:punct:]]", "")))
    gdsc_names <- tolower(gsub(" ","",str_replace_all(df$drug, "[[:punct:]]", "")))
    
    resp <- ccle$resp[[type]]
    cosmics <- unlist(lapply(row.names(resp), function(x) paste0("",as.character(ccle$cl_info$cosmic_id[ccle$cl_info$master_ccl_id == x]))))
    S <- get_data_for_scatter(df = df, results_collected = results_collected, which.cancer.type = type)
    resp <- resp[cosmics %in% row.names(S$RESULTS$DAT$meth),]
    row.names(resp) = cosmics[cosmics %in% row.names(S$RESULTS$DAT$meth)]
    int <- intersect(row.names(S$RESULTS$DAT$meth), row.names(resp))
    resp <- resp[int,]
    meth <- S$RESULTS$DAT$meth[int,]
    dfctrptmp <- list(); for(i in 1:nrow(df)){  
      if((gdsc_names[i] %in% ccle_names) & (df$cancertype[i] == type)){ # if drug is in CCLE
        print(i)
        ccle_drug_id <- ccle$drugs_info$master_cpd_id[ccle_names == gdsc_names[i]]
        tmp_resp <- resp[,as.character(ccle_drug_id)]
        names(tmp_resp) = row.names(resp)
        if(length(which(!is.na(tmp_resp)))>7){ # only at least 10 non-na's
          cpgs <- df$cpg[[i]]
          cpgs <- cpgs[cpgs %in% colnames(meth)]
          dfctrptmp[[i]] <- cpg.assoc_my(beta.val = meth[,cpgs,drop = F], 
                       pheno = tmp_resp/16, 
                       logit.transform = T, 
                       covariates = S$RESULTS$DAT$covariates[int,,drop =F],
                       lrtest = F)
          dfctrptmp[[i]]$index <- i
          dfctrptmp[[i]]$drugctrp <- gdsc_names[i]
          dfctrptmp[[i]]$genectrp <- df$GGene[i]
          dfctrptmp[[i]]$cancertype <- type
        }
      }else{
        next
      }
    }
    dfctrp <- do.call(rbind, dfctrptmp)
    ctrptmp[[type]] <- dfctrp
  }
  ctrp <- do.call(rbind, ctrptmp)
  save(ctrp, file = paste0(savepath, "/ctrp_validation.RData"))
}
