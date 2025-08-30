# Functional module 1: ORA


# -------------------Load data and library -------------

library(plyr)
library(tidyverse)
library(data.table)


dir.create("./Output/M1_ORA", recursive = TRUE, showWarnings = FALSE)


HBM_HMDB_metabolite_ID_fullList <- readRDS("./Data/HBM_HMDB_metabolite_ID_fullList.RDS")
HUBMet_annotation <- readRDS("./Data/HUBMet_annotation.RDS") 
HUBMet_term <- readRDS("./Data/HUBMet_term.RDS")  



#--------------- Functions ----------

## -------- statistical method for ORA and TSA ---------

fisher_enrich <- function(term_met, met_id, n_met_background){
  # Fisher's Exact Test
  term_met <- unique(na.omit(strsplit(term_met, split = ";")[[1]]))   
  overlap <- intersect(term_met, met_id)
  overlap_n <- length(overlap)
  xx <- n_met_background-length(met_id)-length(term_met)+overlap_n
   
  datax <- matrix(c(overlap_n, length(met_id)-overlap_n, 
                    length(term_met)-overlap_n, xx), nrow = 2)
  pp <- fisher.test(datax, alternative = "greater")
  out <- c(overlap_n, paste(overlap, collapse = ";"), pp$p.value,  
           overlap_n/length(met_id) )
  return(out)
}




binomial_enrich <- function(term_met, met_id, n_met_background){
  # Binomial Test
  term_met <- unique(na.omit(strsplit(term_met, split = ";")[[1]]))  
  overlap <- intersect(term_met, met_id)
  overlap_n <- length(overlap)
  pp <- binom.test(x = overlap_n, n = length(met_id), p = length(term_met)/n_met_background, alternative = "greater" )
  out <- c(overlap_n, paste(overlap, collapse = ";"),  pp$p.value,
           overlap_n/length(met_id))
  return(out)
  
}






## --------- background adjustment -------------

background_adjust <- function(dataset, background_hmdb_list){
  
  inner_func <- function(x){
    x <- unique(na.omit(intersect(setdiff(na.omit(unique(strsplit(x, split = ";")[[1]])),c("","NA")),
                                  background_hmdb_list)))
    return(c(paste(x, collapse = ";"), length(x))  )
  }
  
  temp <- as.data.frame(t(as.data.frame(lapply(dataset$HMDB.ID, function(x){inner_func(x)}))))
  colnames(temp) <- c("HMDB.ID", "N_met")
  dataset$HMDB.ID <- as.character(temp$HMDB.ID)
  dataset$N_met <- as.integer(temp$N_met)
  rm(temp)
  return(dataset[which(dataset$N_met >= 5),])
}

 



## ---------Prepare result table------------


table_export <- function(da, database ){
  
  library(tidyverse)
  if(database == "smpdb"){
    colnames(da)[1:3] <- c("SMPDB.Name", "SMPDB.ID", "Pathway.Class")
  }else if(database == "kegg"){
    colnames(da)[1:3] <- c("KEGG.Name", "KEGG.ID", "Pathway.Class")
  } else if(database == "reactome"){
    colnames(da)[1:2] <- c("Reactome.Name","Reactome.ID")
    da <- da[,-c(3,4)]
  }else if(database == "humanGEM"){
    colnames(da)[1] <- "HumanGEM-Subsystem"
    da <- da[,-c(2)]
  }else if(database %in% c("class","subclass","superclass","disease") ){
    colnames(da)[1] <- paste("HMDB",stringr::str_to_title(colnames(da)[1]),sep="-")
  }else if(database %in% c("HUBMet","Drug", "Disease", "Class", "Pathway") ){
    database <- paste("HUBMet", database, sep="_")   
  }
  
  res_output <- da[,which(!(colnames(da)  %in% c("HBM_ID","HMDB.ID","FDR_pri")  ))]
  
  colnames(res_output)[which(colnames(res_output) %in% c("N_met"))] <- "TermSize"
  colnames(res_output)[which(colnames(res_output) %in% c("match_N","hit_n"))] <- "N_Hit"
  colnames(res_output)[which(colnames(res_output) %in% c("match","hit_HBM_ID","hit_HMDB_ID"))] <- "Hit"
  
     
   
  write.table(res_output ,  
              paste("./Output/M1_ORA/",database,"_enrich.txt",sep=""),
              row.names = F, sep = "\t",quote = T)
 
}




##------------------ Main function ----------------

### ----------------- HUBMet database ----------------------
HUBMet_enrich <- function(HBM_list, termclass = "Class" ,method = "fisher", adjp_method = "BH", adjustBackground = "blood", background = NA){
  
  # essential data: HUBMet_term, HUBMet_annotation
  
  # Output: one result table in the folder "Output/M1_ORA/" 
  
  # process HBM ID
  HBM_list <- unique(na.omit(HBM_list[!(HBM_list %in% c("NA",""))]))
  n_input <- length(unique(na.omit(HBM_list)))
  
  res <- HUBMet_term %>%
    filter(TermClass == termclass)
   
  if(adjustBackground %in% c("custom", "blood") ){
    if(adjustBackground == "blood"){
      background <- HUBMet_annotation$HBM_ID
    }
    background <- unique(na.omit(union(HBM_list, 
                                       na.omit(background[!(background %in% c("","NA"))]))))
    dataset <- res
    inner_func <- function(x){
      x <- intersect(setdiff(na.omit(unique(strsplit(x, split = ";")[[1]])),c("","NA")),
                     background)
      return(c(paste(x, collapse = ";"), length(x))  )
    }
    
    temp <- as.data.frame(t(as.data.frame(lapply(dataset$HBM_ID, function(x){inner_func(x)}))))
    colnames(temp) <- c("HBM_ID", "N_met")
    dataset$HBM_ID <- as.character(temp$HBM_ID)
    dataset$N_met <- as.integer(temp$N_met)
    rm(temp)
    res <- dataset[which(dataset$N_met >= 5),]
    
    # Controls the case where there is no entry after background correction N_met >= 5
    if(nrow(res) == 0){
      print("No term has at least 5 metabolites after background adjustment!")
      print("Perform analysis without background adjustment.")
      res <- HUBMet_term
    }
  }
  
  n_met_background <- c(HBM_list, background, strsplit(paste(res$HBM_ID,collapse = ";"),";")[[1]] )
  n_met_background <- length(unique(na.omit(n_met_background[!(n_met_background %in% c("","NA"))])))
  
  if(method == "fisher"){
    res[,c("match_N","match","pvalue","MetaboliteRatio")] <- 
      t(data.frame(lapply(res$HBM_ID, fisher_enrich, met_id = HBM_list, n_met_background)))
  }
  
  if(method == "binomial"){
    res[,c("match_N","match","pvalue","MetaboliteRatio")] <- 
      t(data.frame(lapply(res$HBM_ID, binomial_enrich, met_id = HBM_list, n_met_background)))
  }
  
  res$match_N <- as.integer(res$match_N)
  res <- res[which(res$match_N != 0),]
  
  #  backgroud size 
  n_met_background <- length(unique(strsplit(paste(res$HBM_ID, collapse = ";"), split = ";")[[1]]))
  
  # other parematers
  for(col in c("pvalue","N_met","match_N","MetaboliteRatio" )){
    res[,col] <- as.numeric(res[,col])
  }
  res <- res[order(res$pvalue, decreasing = F),]
  res$BgRatio <- res$N_met/n_met_background
  res$FoldEnrichment <- res$MetaboliteRatio/res$BgRatio
  res$EnrichmentFactor <- (res$match_N/res$N_met)/(length(unique( 
    strsplit(paste(res$match, collapse = ";"),split = ";")[[1]] ))/n_met_background)
  res$p.adjust <- p.adjust(res$pvalue, method = adjp_method  ) 
  
  # export result table
  table_export(res,paste("HUBMet",termclass,sep="_")  )   
  
  return(res)
}




# ------------- Test ------------------

testda <- rio::import("./testda/testData1.txt")

res_ora1 <- HUBMet_enrich(HBM_list = unique(na.omit(testda$metID)),
                         termclass = "Pathway", 
                         adjustBackground = "blood")

res_ora2 <- HUBMet_enrich(HBM_list = unique(na.omit(testda$metID)),
                          termclass = "Class", 
                          adjustBackground = "blood")



# Visualization

source("./Script/visualization.R") # functions for visualization

dotplot_enrich(res_ora1, database = "Pathway")
 
