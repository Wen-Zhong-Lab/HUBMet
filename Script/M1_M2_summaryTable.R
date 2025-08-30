
# Summary table for ORA and MSEA
library(tidyverse)
library(data.table)


sim_matrix_all <- read_rds("./Data/sim_matrix_all.RDS")  

 
# chosen_database is a list:
# options: "HUBMet_Pathway","HUBMet_Drug","HUBMet_Disease","HUBMet_Class","class","subclass","superclass","disease","humanGEM","kegg","reactome","smpdb"

# ----------------- summary table for ORA ------------------------

fuzzy_cluster_m1 <- function(chosen_database){
  
  results_list <- list()
 
  for (i in chosen_database) {
    
    e1 <- fread(paste0("./Output/M1_ORA/", i, "_enrich.txt"), sep = "\t", header = TRUE, fill = TRUE)
    if (i == "smpdb"){
      colnames(e1)[3] <- "Class"
      e1$Class <- paste("smpdb_",e1$Class,sep="")
      colnames(e1)[2] <- "Supplement"
      e1$Supplement <- paste("SMPDB_ID:",e1$Supplement,sep="")
    }
    else if (i == "kegg"){
      colnames(e1)[3] <- "Class"
      e1$Class <- paste("kegg_",e1$Class,sep="")
      colnames(e1)[2] <- "Supplement"
      e1$Supplement <- paste("KEGG_ID:",e1$Supplement,sep="")
    }
    else if (i == "reactome"){
      colnames(e1)[2] <- "Supplement"
      e1$Supplement <- paste("Reactome_id:",e1$Supplement,sep="")
      e1$Class <- i
    }
    else if (i %in% c("HUBMet_Pathway","HUBMet_Drug","HUBMet_Disease","HUBMet_Class")){
      e1$Term <- sapply(strsplit(as.character(e1$Term), ":"), `[`, 1)
      colnames(e1)[2] <- "Class"
      e1$Class <- paste("HUBMet_",e1$Class,sep="")
      e1$Supplement <- NA
    }
    else {
      colnames(e1)[1] <- "Name"
      e1$Class <- i
    }
    colnames(e1)[1] <- "Name"
    e1$Name <- paste(e1$Name, e1$Class, sep = "_") 
    results_list[[i]] <- e1  
  }
  
  results <- rbindlist(results_list, fill = TRUE) %>% as.data.frame()  
  
  results <- results[order(results$p.adjust), ]
  results$order <- 1:nrow(results)
  results$group <- 1:nrow(results)
  
  if (nrow(results) > 1){
    for (i in results$order[1:(nrow(results)-1)]) {
      if (results$order[i]==results$group[i]){
        for (j in (i+1):nrow(results)) {
          if (results$order[j]==results$group[j]) {
            if (sim_matrix_all[results$Name[i],results$Name[j]] > 0.3) {
              results$group[j] <- i
            }
          }
        }
      }
    }
    results$main <- results$order==results$group
    grp_idn <- sim_matrix_all[results$Name[results$main==TRUE],results$Name[results$main==FALSE]]
    for (i in colnames(grp_idn)) {
      more_sig <- results$Name[results$p.adjust <= results$p.adjust[results$Name==i] & results$main==TRUE]
      groupnumber <- results[results$Name == more_sig[which.max(grp_idn[more_sig,i])],"group"]
      results[results$Name == i,"group"] <- unique(groupnumber)
    }
    results <- results[order(results$group), ]
    results$Group <- as.integer(as.factor(results$group))
    results <- results[,c(1,16,2:12,15)]
    colnames(results)[c(1,3)] <- c("Term.Name","TermClass")
    results[which(is.na(results$Supplement)),"Supplement"] <- ""
    write.table(results ,    "./Output/M1_ORA/summary_table_enrich.txt",row.names = F, sep="\t", quote = T)
    return(results)
  } else {
    results$main <- results$order==results$group
    results$Group <- as.integer(as.factor(results$group))
 
    results <- results[,c(1,16,2:12,15)]
    colnames(results)[c(1,3)] <- c("Term.Name","TermClass")
    results[which(is.na(results$Supplement)),"Supplement"] <- ""
    write.table(results[,!colnames(results) %in% c("Hit_link","Hit_color")],  
                "./Output/M1_ORA/summary_table_enrich.txt",row.names = F, sep="\t", quote = T)
    return(results)
  }
}


# ----------------- summary table for MSEA ------------------------

fuzzy_cluster_m2 <- function(chosen_database){
  
  results_list <- list() 
  for (i in chosen_database) {
    e1 <- read.delim(file=paste("./Output/M2_MSEA/",i,"_MSEA.txt",sep=""),sep="\t",header=TRUE,check.names=FALSE,fill = TRUE)
    if (i == "smpdb"){
      colnames(e1)[3] <- "Class"
      e1$Class <- paste("smpdb_",e1$Class,sep="")
      colnames(e1)[2] <- "Supplement"
      e1$Supplement <- paste("SMPDB_ID:",e1$Supplement,sep="")
    }
    else if (i == "kegg"){
      colnames(e1)[3] <- "Class"
      e1$Class <- paste("kegg_",e1$Class,sep="")
      colnames(e1)[2] <- "Supplement"
      e1$Supplement <- paste("KEGG_ID:",e1$Supplement,sep="")
    }
    else if (i == "reactome"){
      colnames(e1)[2] <- "Supplement"
      e1$Supplement <- paste("Reactome_id:",e1$Supplement,sep="")
      e1$Class <- i
    }
    else if (i %in% c("HUBMet_Pathway","HUBMet_Drug","HUBMet_Disease","HUBMet_Class")){
      e1$Term <- sapply(strsplit(as.character(e1$Term), ":"), `[`, 1)
      colnames(e1)[2] <- "Class"
      e1$Class <- paste("HUBMet_",e1$Class,sep="")
      e1$Supplement <- NA
    }
    else {
      e1$Class <- i
    }
    colnames(e1)[1] <- "Name"
    e1$Name <- paste(e1$Name, e1$Class, sep = "_")
    results_list[[i]] <- e1  
  }
  results <- rbindlist(results_list, fill = TRUE) %>% as.data.frame()  
  
  results <- results[order(results$FDR), ]
  results$order <- 1:nrow(results)
  results$group <- 1:nrow(results)
  if (nrow(results) > 1){
    for (i in results$order[1:(nrow(results)-1)]) {
      if (results$order[i]==results$group[i]){
        for (j in (i+1):nrow(results)) {
          if (results$order[j]==results$group[j]) {
            if (sim_matrix_all[results$Name[i],results$Name[j]] > 0.3) {
              results$group[j] <- i
            }
          }
        }
      }
    }
    results$main <- results$order==results$group
    grp_idn <- sim_matrix_all[results$Name[results$main==TRUE],results$Name[results$main==FALSE]]
    for (i in colnames(grp_idn)) {
      more_sig <- results$Name[results$FDR <= results$FDR[results$Name==i] & results$main==TRUE]
      groupnumber <- results[results$Name == more_sig[which.max(grp_idn[more_sig,i])],"group"]
      results[results$Name == i,"group"] <- unique(groupnumber)
    }
    results <- results[order(results$group), ]
    results$Group <- as.integer(as.factor(results$group))
    
     
    results <- results[,c(1,14,2:10,13)] 
    colnames(results)[c(1,3)] <- c("Term.Name","TermClass")
    results[which(is.na(results$Supplement)),"Supplement"] <- ""
    write.table(results ,  "./Output/M2_MSEA/summary_table_MSEA.txt",row.names = F, sep="\t",quote = T) 
    return(results)
  } else {
    results$main <- results$order==results$group
    results$Group <- as.integer(as.factor(results$group))
     
    results <- results[,c(1,14,2:10,13)]  
    results$Name <- sapply(strsplit(as.character(results$Name), "_"), `[`, 1)
    colnames(results)[c(1,3)] <- c("Term.Name","TermClass")
    results[which(is.na(results$Supplement)),"Supplement"] <- ""
    write.table(results ,   
                "./Output/M2_MSEA/summary_table_MSEA.txt",row.names = F, sep="\t",quote = T)
     return(results)
  }
}


 
# ----------------Test --------------------

summary <- fuzzy_cluster_m1(c("HUBMet_Pathway","HUBMet_Class" ))

summary <- fuzzy_cluster_m2(c("HUBMet_Pathway","HUBMet_Drug","HUBMet_Disease","HUBMet_Class" ))
