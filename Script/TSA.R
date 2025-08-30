# Functional module 3 TSA
# -------------------Load data and library -------------

library(plyr)
library(tidyverse)

HUBMet_annotation <- readRDS("./Data/HUBMet_annotation.RDS")    
HUBMet_tissue_criteria <- readRDS("./Data/HUBMet_tissue_criteria.RDS")  
HUBMet_tissue_term <- readRDS("./Data/HUBMet_tissue_term.RDS")  

 
dir.create("./Output/M3_TSA", recursive = TRUE, showWarnings = FALSE)

# -------------------Function --------------------

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


## ------------- Main function ---------------

HUBMet_Tissue <- function(HBM_list, method = "fisher", adjp_method = "BH", adjustBackground = "blood", background = NA,
                          job_id="job_id"){
  library(ggrepel)
  
  # essential data : HUBMet_tissue_term, HUBMet_annotation, HUBMet_tissue_criteria
  # output: a result table, one pie plot and one bar plot in the folder 'Output/M3_TSA'
  
  # process HBM ID
  HBM_list <- unique(na.omit(HBM_list[!(HBM_list %in% c("NA",""))]))
  n_input <- length(unique(na.omit(HBM_list)))
  
  # if no metabolite has been identified in HBMD, print the sentence and end the function
  if(length(HBM_list) == 0){
    return("No metabolite has been mapped to HBMD!")
  }
  
  # if no metabolite has tissue specificity, to avoid print sentences during background adjustment
  temp <- HUBMet_tissue_criteria[which(HUBMet_tissue_criteria$HBM_ID %in% HBM_list),]
  temp <- setdiff(temp$Reliability,c("Low", "Low Specificity", "Unknown"))
  if(length(temp) == 0){
    adjustBackground == "no"
  }
  
  # start Tissue enrichment analysis
  res <- HUBMet_tissue_term[which(HUBMet_tissue_term$N_met >= 3),]
  
  # background adjustment 
  if(adjustBackground %in% c("custom", "blood","all") ){
    if(adjustBackground != "custom"){
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
    res <- dataset[which(dataset$N_met >= 3),]
    
    if(nrow(res) == 0){
      print("No term has at least 3 metabolites after background adjustment!")
      print("Perform analysis without background adjustment.")
      res <- HUBMet_tissue_term[which(HUBMet_tissue_term$N_met >= 3),]
    }
  }
  
  n_met_background <- c(HBM_list, background, strsplit(paste(res$HBM_ID,collapse = ";"),";")[[1]] )
  n_met_background <- length(unique(na.omit(n_met_background[!(n_met_background %in% c("","NA"))])))
  
  if(method == "fisher"){
    res[,c("match_N","match","pvalue","MetaboliteRatio")] <- 
      t(data.frame(lapply(res$HBM_ID, fisher_enrich, met_id = HBM_list, n_met_background)))
  }else if(method == "binomial"){
    res[,c("match_N","match","pvalue","MetaboliteRatio")] <- 
      t(data.frame(lapply(res$HBM_ID, binomial_enrich, met_id = HBM_list, n_met_background)))
  }
  
   
  res <- as.data.frame(res) %>%  
    mutate( match_N = as.integer( match_N)) %>% 
    filter(match_N != 0)
  
  n_met_background <- length(unique(strsplit(paste(res$HBM_ID, collapse = ";"), split = ";")[[1]]))
  
  for(col in c("pvalue","N_met","match_N","MetaboliteRatio" )){
    res[,col] <- as.numeric(res[,col])
  }
  
  res <- res %>%
    arrange(pvalue) %>%
    mutate(
      BgRatio = N_met / n_met_background,
      FoldEnrichment = MetaboliteRatio / BgRatio,
      EnrichmentFactor = (match_N / N_met) / 
        (length(unique(unlist(strsplit(paste(match, collapse = ";"), ";")))) / n_met_background),
      p.adjust = p.adjust(pvalue, method = adjp_method)
    )
  
  res_output <- res[,which(!(colnames(res)  %in% c("HBM_ID","HMDB.ID","FDR_pri")  ))]
  
  colnames(res_output)[which(colnames(res_output) %in% c("N_met"))] <- "TermSize"
  colnames(res_output)[which(colnames(res_output) %in% c("match_N","hit_n"))] <- "N_Hit"
  colnames(res_output)[which(colnames(res_output) %in% c("match","hit_HBM_ID","hit_HMDB_ID"))] <- "Hit"
 
  fn = paste(c("./Output/M3_TSA/",job_id,"_HUBMet_tissue_specificity.txt"),collapse ="")
  write.table(res_output ,         
              fn,row.names = F, sep = "\t",quote = T)
  
  
  # visualization
  
  da <- res[which(res$p.adjust < 0.05),!(colnames(res) %in% c("BTO_ID"))]  
  
  if(nrow(res) == 0){
    tissue_bar <- ggplot()+theme(panel.background = element_rect(fill='transparent'))
    print("No metabolite has been identified with tissue specificity!")
  }else if(nrow(da) == 0){
    tissue_bar <- ggplot()+theme(panel.background = element_rect(fill='transparent'))
    file.create( paste0( "./Output/M3_TSA/",job_id,"_tissue_bar_pie_2_empty.csv")  ) 
  }else{
    
    da <- da[,c(1,4,6,11)]  
    da[which(da$p.adjust < 0.05),"sig"] <- "*"
    da[which(da$p.adjust < 0.01),"sig"] <- "**"
    da[which(da$p.adjust < 0.001),"sig"] <- "***"
    da[which(da$p.adjust < 0.0001),"sig"] <- "****"
    da[is.na(da$sig),"sig"] <- ""
    
    da <- da[order(-da$match_N),] %>%
      mutate(Tissue = factor(Tissue, levels = Tissue),
             label = paste(sig, match_N, sep="\n"))  
    
    
    # control the range of y-axis
    if(max(da$match_N) > 200){
      yadd <- ceiling(da$match_N*0.2)  
    }else if(max(da$match_N) > 100){
      yadd <- 15+5   
    }else if(max(da$match_N) < 50){
      yadd <- 7+5   
    }else{
      yadd <- 10+5  
    }
    
    # control the font size on bar
    if(nrow(da) > 25){
      xsize <- 9
      nsize <- 2.5 
    }else{
      xsize <- 12
      nsize <- 3.5 
    }
    
    tissue_bar <- ggplot(da,aes(x=Tissue, y=match_N))+
      geom_bar(stat = "identity",aes(alpha=match_N), fill="#8B668B" )+ 
      geom_text(aes(label=label,x= Tissue ,y=match_N), vjust = 0, hjust = 0.5,size = nsize) +  
      scale_y_continuous(expand = c(0,0),limits = c(0,max(da$match_N+yadd)))+
      scale_alpha_continuous(range = c(0.45,1))+
      labs(x=expression("*" * italic("adj.P") * " < 0.05; **" * italic("adj.P") * " < 0.01; *** " * italic("adj.P") * " < 0.001; **** " * italic("adj.P") * " < 0.0001"),
           y="Number of metabolite")+
      theme_classic() +
      theme(legend.position = "none",
            axis.line = element_line(color="black", linewidth = 0.5),
            axis.title.y  =element_text(size=12, colour = "black"), 
            axis.title.x =element_text(size=10, colour = "black"),
            axis.text.x =element_text(size=xsize, colour = "black",angle = 30,hjust = 1), 
            axis.text.y = element_text(size=12, colour = "black"),
            plot.margin = margin(rep(2.5,4),  unit = "line"))
    
 
    
  }
  
  
  
  # tissue specificity 2
  
  state_color <- c("#8B8970", "#CDC8B1", "#EEE8CD")
  names(state_color) <- c("Tissue Enriched", "Low Specificity", "Unknown")
  
  stat_color_df <- as.data.frame(state_color) %>% 
    rownames_to_column(var="state")
  
  
  da <- HUBMet_annotation[which(HUBMet_annotation$HBM_ID %in% HBM_list),]
  da$state <- da$Reliability.Tissue
  da[which(!da$Reliability.Tissue %in% c("Unknown","Low Specificity")),"state"] <- "Tissue Enriched"
  
  
  
  a2 <- da[,c("state","HBM_ID","Name")]  
  
  a2 <- table(da$state) %>% as.data.frame() %>% 
    mutate( state = factor(Var1, levels = c("Tissue Enriched", "Low Specificity", "Unknown")),
            label = paste0(state," (n=",Freq, ")"),
            text = paste0(state,"\n(n=",Freq, ")"),
            perc  = sprintf("%.2f%%", Freq / sum(Freq) * 100),
            text = paste0(label,"\n",perc) )
  
   
  pie_tissue_state <- ggplot(a2, aes(x="", y=Freq, fill=state)) +
    geom_bar(stat="identity", width=0.25, color="white") +
    coord_polar("y", start=0)+
    geom_text(aes(label = text ), size=3.5,
              position = position_stack(vjust = 0.5)) + 
    scale_fill_manual(values = state_color)+
    theme_classic()+
    theme(plot.background = element_rect(fill='transparent', color=NA),
          panel.background = element_rect(fill='transparent'),
          legend.position = "none",
          axis.title.x =element_blank(),
          axis.title.y =element_blank(),
          axis.text.x =element_blank(),
          axis.text.y =element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank())
  
  
  
  
  #  layout
  p_tissue <- cowplot::plot_grid(pie_tissue_state,tissue_bar,rel_widths = c(0.5,1.25))   
  
  if(nrow(res) == 0 | nrow(res[which(res$p.adjust < 0.05),])==0){
    p_tissue <- pie_tissue_state
  }
  
  fn = paste(c("./Output/M3_TSA/",job_id,"_tissue_bar_pie.pdf"),collapse ="")
  ggsave(fn,p_tissue, width = 8, height = 4.5)  
  
  
}


# ---------- Test -------------

testda <- rio::import("./testda/testData1.txt")

res <- HUBMet_Tissue(na.omit(unique(testda$metID)),adjustBackground = "blood",job_id = "test42")


