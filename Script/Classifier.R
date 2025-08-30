# Module classifier

# ------------- Load data and library ------------

library(tidyverse)

 
dir.create("./Output/M0_Classifier", recursive = TRUE, showWarnings = FALSE)


HUBMet_annotation <- readRDS("./Data/HUBMet_annotation.RDS")  


# -------------------- Main function -----------------------

HUBMet_basic <- function(HBM_list,job_id="job_id"){
  
  # input HBM_list is the list of met_ID in the mapping table
  # output in the folder '/Output/Classifier': result table of annotation ,  bar plot of category
  
  
  # essential data:   HUBMet_annotation
  library(ggplot2)
  library(cowplot)
  
  HBM_list <- unique(na.omit(HBM_list[!(HBM_list %in% c("NA",""))]))
  temp <- HUBMet_annotation[which(HUBMet_annotation$HBM_ID %in% HBM_list),]
  
  if(nrow(temp) == 0){
    return("No metabolite is mapped to HUBMet!")
  }else{
    print(paste(nrow(temp)," metabolite(s) are mapped to HUBMet!",sep=""))
  }
  
  
  temp_need <- temp  
  temp_need[is.na(temp_need)] <- ""
  
  
  colnames(temp_need)[c(1,3,4,5)] <-  
    c("HUBMet ID","Category", "Tissue.Relevance","Tissue.Relevance.Reliability") 
  # export table
  write.table(temp_need,
              paste("./Output/M0_Classifier/",job_id,"_basic_annotation.txt",sep = ""),
              sep="\t", row.names = F, quote = F)
  rm(temp_need)
  
  
  # category
  color9_new <- c("#ffa510", "#45aee8", "#EB99C2", "#f6c9aa", "#f57070", 
                  "#add547", "#41b7ac", "#ad85d5", "#C3DEE0")
  
  Category9_level <- c("Lipid","Xenobiotics","Amino Acid","Peptide","Carbohydrate",
                       "Nucleotide","Cofactors and Vitamins","Energy","Others")
  names(color9_new) <- Category9_level
  
  
  color9_new_df <-  data.frame("color"=color9_new,"category"=Category9_level)
  
  
  da2 <- temp %>% 
    dplyr::select(Category, Name, HBM_ID) %>% 
    mutate( id= paste(  "https://hubmet.app.bio-it.tech/metabolite/",HBM_ID,sep=""),
            name = gsub(";"," ",Name)) %>% 
    group_by(Category ) %>% 
    summarise(links=paste(id,collapse = ";"),
              labels=paste(name,collapse = ";"),
              n=length(unique(HBM_ID)),
              .groups = "drop") %>%
    mutate( term=Category,
            label = paste0(Category, " (n=", n, ")"),
            Category = factor(Category, levels = Category),
            perc  = sprintf("%.2f%%", n / sum(n) * 100),
            text = paste0(label, "\n", perc)) %>%
    arrange(desc(n)) %>% 
    mutate(Category = factor(Category, levels = Category))
  
  
  pie_Category <- ggplot(da2, aes(x="", y=n, fill=Category)) +
    geom_bar(stat="identity", width=1, color="white", linewidth=0.5) +
    coord_polar("y", start=0)+
    scale_fill_manual(values = color9_new, breaks=da2$Category, labels = da2$label,
                      name="Category")+  
    theme_classic()+
    theme(legend.position = "right",
          axis.title.x =element_blank(), 
          axis.title.y =element_blank(),
          axis.text.x =element_blank(), 
          axis.text.y =element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          legend.text=element_text(size=11, colour = "black"), 
          legend.title=element_text(size=11, colour = "black"))
  
  # export figures
  ggsave(paste("./Output/M0_Classifier/",job_id,"_mainclass_pie.pdf",sep = ""),
         pie_Category, width = 7, height = 3)

}


# -------------------- Test  ----------------------

set.seed(42)

x <- HUBMet_annotation$HBM_ID[sample(1:3950,500,replace = F)]

HUBMet_basic(HBM_list = x, job_id = "test42")

