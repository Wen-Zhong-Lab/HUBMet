# Visualization for ORA and MSEA

library(plyr)
library(tidyverse)
library(pathview)
library(SBGNview)
library(SBGNview.data)
library(ggrepel)
library(scales)

HUBMet_annotation <- readRDS("./Data/HUBMet_annotation_250827.RDS") 
HMDB_ID_unique_name <- read_rds("./Data/HMDB_ID_unique_name_1028.RDS")  
HUBMet_term <- readRDS("./Data/HUBMet_term_1028.RDS") 
# --------------- ORA --------------
## ----------Dot plot in ORA-------------------------------


# 2507 dotplot_enrich 主要修改： 输出csv

dotplot_enrich <- function(da, database, 
                           smpdb_class = "Metabolic", 
                           kegg_class = "Amino acid metabolism" ){  
  
  da <- da %>%
    mutate(pvalue = as.numeric(pvalue),
           p.adjust = as.numeric(p.adjust)) %>%
    filter(p.adjust < 0.05)
  
  
  if(nrow(da[which(!is.na(da$pvalue) & da$p.adjust < 0.05),]) == 0){  # change 1206
    p <- "No significantly enriched term (adjusted p-value < 0.05) in this section!"  # change 1206
    return(p)
  }
  
  theme_dotplot <- theme_classic(base_size = 13) +  # 1028 change theme_bw() to theme_classic
    theme(plot.title = element_text(hjust = 0.5,size=rel(1), colour = "black", face = "bold"),
          legend.position = "right",
          axis.line = element_line(linewidth = 0.5),
          legend.text = element_text(size=rel(0.85), colour = "black"),
          legend.title = element_text(size=rel(0.9), colour = "black"),
          axis.title =element_text(size=rel(1), colour = "black"), 
          axis.text =element_text(size=rel(1), colour = "black") )
  color_legend <- "p.adjust"
  
  if(database == "smpdb"){
    da <- da[which(da$Pathway.Class == smpdb_class),]
    ti = paste(c("SMPDB pathway", " (", smpdb_class,")"),collapse = "")
    fn = paste("./Output/M1_ORA/",database,"_",smpdb_class,"_dotplot.pdf",sep = "")
  }else if(database == "kegg"){
    da <- da[which(da$Pathway_Class ==  kegg_class),]
    ti =  paste(c("KEGG pathway", " (", kegg_class,")"),collapse = "")
    fn = paste("./Output/M1_ORA/",database,"_",kegg_class,"_dotplot.pdf",sep = "")
  }else if(database  %in% c("Class","Pathway","Disease","Drug")){  # 修改，将HUBMet四类term分别作为选项进行计算
    da <- da[which(da$TermClass ==  database),]
    ti = paste(c(database,"HUBMet" ),collapse = "-")
    fn = paste("./Output/M1_ORA/", "HUBMet_",database,"_dotplot.pdf",sep = "")
  }else{
    fn = paste("./Output/M1_ORA/",database,"_dotplot.pdf",sep = "")
    if(database == "humanGEM"){
      ti = "HumanGEM-Subsystem" # 1028
    }else if(database == 'reactome'){
      ti = "Reactome pathway"
    }else if(database %in%  c("disease","class","subclass","superclass")){
      ti = paste("HMDB",stringr::str_to_title(colnames(da)[1]),   sep = "-") # 1028
    }   
  } 
  
  
  da <- da[1:20,] %>% filter(!is.na(BgRatio)) 
  
  
  if(nrow(da) == 0){
    p <- "No enriched term or adjusted p-value < 1 in this section!"
  }else{
    danrow <- nrow(da) # 用于控制输出的pdf的height
    colnames(da)[1] <- "term"
    
    # 创建一个颜色渐变映射函数（和 ggplot 的 low-high 对应）
    color_fun <- scales::col_numeric(palette = c("#e04c4c", "#fdd819"), domain = range(da$p.adjust, na.rm = TRUE))
    
    
    da <- da %>%
      mutate(
        pvalue = as.numeric(pvalue),
        p.adjust = as.numeric(p.adjust),
        match_N = as.integer(match_N),
        term = factor(term, levels = term[order(EnrichmentFactor )]),
        color = if (min(p.adjust) == max(p.adjust) & min(p.adjust) > 0.05) {
          "#BEBEBE"
        } else {
          color_fun(p.adjust)
        })
    
    
    
    p <- ggplot(da,aes(x=EnrichmentFactor,y=term)) + 
      geom_point(mapping = aes(size = match_N, color = p.adjust))+
      scale_y_discrete(labels=function(x) str_wrap(x, width=32,whitespace_only = FALSE))+
      guides(color = guide_colorbar(order = 2 ), size = guide_legend(order = 1))+
      labs(x="Enrichment Factor", title=element_blank(), 
           size = "Count", color= "p.adjust",  y = ti)+
      theme_dotplot
    # 基于校正后p值范围修改颜色
    if(min(da$p.adjust) == max(da$p.adjust) & min(da$p.adjust) > 0.05){
      p <- p + scale_color_gradient(high="#BEBEBE", low="#BEBEBE")
    }else{
      # change 1206
      p <- p + scale_color_gradient(  low="#e04c4c",high ="#fdd819" )
    }
    
    # 基于绘制的条目数量调整pdf的高度height, 宽度with保持为7  
    if(danrow < 15){
      ggsave(filename = fn,p,width = 7,height = 5)
    }else{
      ggsave(filename = fn,p,width = 7,height = 7/20*danrow)  
    }
    
      
  }
  return(p)
}



# ----------------- MSEA ---------------
 
## ----------------------- Data preparation for MSEA -------------------
solve_redundant <- function(hmdb_list, hmdb_list_value, choice = "mean"){
  da <- data.frame("hmdb_list" = hmdb_list,"hmdb_list_value"= hmdb_list_value)
  if(length(unique(hmdb_list)) < length(hmdb_list)){
    a <- data.frame(table(hmdb_list))
    a <- a[which(a$Freq>1),]
    a$hmdb_list  <- as.character(a$hmdb_list)
    if(choice == "mean"){
      for( i in a$hmdb_list){
        da[which(da$hmdb_list == i),"hmdb_list_value"] <- mean(da[which(da$hmdb_list == i),"hmdb_list_value"])
      }
      
    }else if(choice == "max"){
      for( i in a$hmdb_list){
        da[which(da$hmdb_list == i),"hmdb_list_value"] <- max(da[which(da$hmdb_list == i),"hmdb_list_value"])
      }
      
    }else if(choice == "min"){
      for( i in a$hmdb_list){
        da[which(da$hmdb_list == i),"hmdb_list_value"] <- min(da[which(da$hmdb_list == i),"hmdb_list_value"])
      }
      
    }
    da <- unique(da)
    
  }
  return(da)
}





## --------------Bar plot in MSEA-------------------------

# barplot_msea 主要修改： export csv

barplot_msea <- function(res_msea, database){
  library(tidyverse)
  top_N=10
  # 处理res_msea无数据的情况，在结果页面提示下列return的句子
  da <- res_msea %>% # change 1206
    mutate(FDR = as.numeric(FDR),
           FDR_pri = as.numeric(FDR_pri),
           NES = as.numeric( NES)) %>% # change 1206
    filter(FDR < 0.25)   # change 1206
  
  if(nrow(da) == 0){
    return("No significantly enriched metabolite set (FDR < 0.25)!")  # change 1206
  }
  
  rm(res_msea)
  
  theme_bar_msea <- theme_classic(base_size = 13) +
    theme(plot.title = element_text(hjust = 0.5,size=rel(1), colour = "black", face = "bold"),
          legend.position = "right",
          legend.text = element_text(size=rel(0.85), colour = "black"),
          legend.title = element_text(size=rel(0.9), colour = "black"),
          axis.line.y = element_blank(),  
          axis.line.x = element_line(linewidth = 0.5,color="black"),
          axis.ticks.y = element_blank(),
          axis.title  =element_text(size=rel(1), colour = "black"),  
          axis.text.x =element_text(size=rel(1), colour = "black"), 
          axis.text.y = element_blank())
  
  fn=database
  
  if(database == "humanGEM"){
    ti = "HumanGEM-Subsystem" # 1028
  }else if(database == 'reactome'){
    ti = "Reactome pathway"
  }else if(database == 'disease'){
    ti = "HMDB-Disease"
  }else if(database %in% c("Class","Pathway","Disease","Drug") ){  # 1028
    ti = paste(c(database,"HUBMet" ),collapse = "-")
    fn = paste("HUBMet",database, sep="_")
  }else if(database == "smpdb"){
    ti = "SMPDB pathway"
  }else if(database == "kegg"){
    ti= "KEGG pathway"
  }
  
  
  colnames(da)[1] <- "term"
  da <- da[order(da$FDR_pri,-abs(da$NES)),] %>% # change
    mutate(logp_FDR = -log10(FDR) ,
           hit_n = as.integer(hit_n)) %>%
    slice_head(n = top_N) %>%
    filter(!is.na(FDR))  
  
  
  if(nrow(da) == 0){
    p = "No enriched metabolite set!"
  }else{
    danrow <- nrow(da) # a control pdf file height   
    wrap_width <- 75
    if( max(da$NES)*min(da$NES) < 0 ){
      comp <- max(c(abs(max(da$NES)), abs(min(da$NES))))
      x_lim <- c(-ceiling(comp), ceiling(comp))
      wrap_width <- 40 
    }else if( max(da$NES) > 0){
      x_lim <- c(0, max(da$NES))
    }else{
      x_lim <- c(min(da$NES), 0 )
    }
    
    da <- da %>%
      mutate(termchr=term,
             term = str_wrap(term, width = wrap_width)) %>% # 移动位置在control length之后
      arrange(desc(NES)) %>%
      mutate( term = factor(term, levels = rev(term)))
    
    
    if( min(da$NES) <= 0 & max(da$NES) <= 0 ){
      da$term <- factor(da$term, levels = da$term[1:length(da$term)]) 
    }
    
    
    p <- ggplot(da,aes(x=NES,y=term, fill = FDR)) + 
      geom_bar(stat="identity",position = "dodge") +
      geom_text(subset(da, NES>=0), mapping = aes(label = term, x=0.01*max(NES)),
                hjust = 0, size = 4, color = "black", alpha = 1)+
      geom_text(subset(da, NES<0), mapping = aes(label = term, x=0.01*min(NES)),
                hjust =1, size = 4, color = "black", alpha = 1)+
      geom_vline(aes(xintercept=0),linewidth=0.5,col="black")+ # change linetype, specify linewidth
      scale_y_discrete(labels=function(x) str_wrap(x, width=30))+
      scale_x_continuous(expand = c(0,0), limits = x_lim)+ # add
      labs(x="Normalized enrichment score", title="Metabolite set enrichment analysis (MSEA)", 
           color="FDR",size = "Hit count", y = ti)+
      theme_bar_msea
    
    if(min(da$FDR) == max(da$FDR) & min(da$FDR) > 0.25){
      p <- p + scale_fill_gradient(high="grey", low="grey")
    }else{ 
      p <- p + scale_fill_gradient(  low="#e04c4c",high ="#fdd819" ) # change 1206
    }
    
    
    # export figure 
    ggsave(filename = paste("./Output/M2_MSEA/",fn,"_barplot.pdf",sep = ""),p,width = 7,height = 5)  # 0404 add
    
    
    # 2507 csv
    
    color_fun <- scales::col_numeric(palette = c("#e04c4c", "#fdd819"), domain = range(da$FDR, na.rm = TRUE))
    da <- da %>%
      mutate( 
        color = if (min(FDR) == max(FDR) & min(FDR) > 0.25) {
          "#BEBEBE"
        } else {
          color_fun(FDR)
        })  
    
     
  }
  #return(p)
}





## ----------------- Heatmap in MSEA --------------------------

msea_heatmap <- function(res_msea, hmdb_list, hmdb_list_value, database){
  top_N = 10
  choice = "mean"
  library(tidyverse)
  # 处理res_msea中无数据的情况
  da <- res_msea %>% # change 1206
    mutate(FDR = as.numeric(FDR),
           FDR_pri = as.numeric(FDR_pri)) %>%  
    filter(FDR < 0.25) 
  
  if(nrow(da) == 0){
    return("No significantly enriched metabolite set (FDR < 0.25)!")
  }
  
  rm(res_msea) 
  da <- da[order(da$FDR_pri, decreasing = F),]
  if(nrow(da)>= top_N){
    da <- da[1:top_N,]
  }else{
    top_N = nrow(da)
  }
  
  # 基于数据库确定pathway名称
  if(database == "kegg"){
    path_name_index <- "KEGG_pathway_name"
    path_id_index <- "KEGG_pathway_ID"
  }else if(database == "smpdb"){
    path_name_index <- "Pathway.Name"
    path_id_index <- "SMP.ID"
  }else if(database == "humanGEM"){
    path_name_index <- "SUBSYSTEM"
    path_id_index <- "SUBSYSTEM"
  }else if(database == "disease"){
    path_name_index <- "disease"
    path_id_index <- "disease"
  }else if(database == "reactome"){
    path_name_index <- "reactome_name"
    path_id_index <- "reactome_id"
  }else if(database %in% c("Class","Pathway","Disease","Drug")){  # HBM的适应化处理
    path_name_index <- "Term"
    path_id_index <- "Term"
    colnames(da)[c(3,6)] <- c("HMDB.ID", "hit_HMDB_ID")
  }
  
  
  hmdb_need <- data.frame("set" = NA, "set_name" = NA, "hmdb" = NA)
  for(i in 1:top_N){
    list_met <- strsplit(da$hit_HMDB_ID[i], split = ";")[[1]]
    temp <- data.frame("set" = rep(da[i,path_id_index], times = length(list_met)),
                       "set_name" = rep(da[i,path_name_index], times = length(list_met)),
                       "hmdb"=list_met)
    hmdb_need <- rbind(hmdb_need, temp )
  }
  hmdb_need <- unique(hmdb_need[which(!is.na(hmdb_need$set)),])
  met_da <- solve_redundant(hmdb_list, hmdb_list_value, choice) %>% 
    filter(hmdb_list %in% unique(hmdb_need$hmdb)) %>% 
    mutate(hmdb_list_value_abs = abs(as.numeric(hmdb_list_value))) %>%
    arrange(desc(hmdb_list_value_abs))
  
  
  max_bk <- ceiling(abs(met_da$hmdb_list_value[1]))  # change 0409
  
  # 修改，控制可视化横轴代谢物数量  
  if(nrow(met_da) > 20){
    met_da <- met_da[1:20,1:2]
  }else{
    met_da <- met_da[,1:2]
  }
  
  # 修改为代谢物名称
  # HBM的适应化处理
  if(database %in% c("Class","Pathway","Disease","Drug")){
    met_da <- left_join(met_da, HUBMet_annotation[,c("HBM_ID","Name")], 
                        by =c("hmdb_list"="HBM_ID") )
    colnames(met_da)[3] <- "name"
  }else{
    met_da <- left_join(met_da, HMDB_ID_unique_name[,c("HMDB_ID_unique","name")], 
                        by =c("hmdb_list"="HMDB_ID_unique") )
  }
  
  
  da_plot <- dplyr::left_join(met_da, hmdb_need , by = c("hmdb_list"="hmdb")) %>% 
    filter(!is.na(hmdb_list)) 
  if(database %in% c("humanGEM","disease","Class","Pathway","Disease","Drug")){ # 1103 change
    da_plot$pathway_name <- as.character(da_plot$set_name)
  }else{
    da_plot$pathway_name <- as.character(unlist(apply(da_plot[,c("set","set_name")],1, 
                                                      function(x){paste(c(x[2]," (", x[1], ")"), collapse = "")})))
  }
  
  da_plot2 <- da_plot %>%
    distinct() %>%
    dplyr::select(pathway_name, name, hmdb_list_value) %>%
    pivot_wider(names_from = name, values_from = hmdb_list_value, values_fill = 0) %>% 
    as.data.frame()
  
  rownames(da_plot2) <- da_plot2$pathway_name
  
  da_plot2 <- da_plot2[,-c(1),drop=FALSE]  # 添加drop == FALSE, 防止只剩一列时被转换成列表
  
  
  # prepare for the plot
  # 不控制代谢物名称长度, 如果设置名称长度限制，横轴代谢物名称会重叠 
  rownames(da_plot2) <- str_wrap(rownames(da_plot2), width = 40)
  
  bk <- unique(c(seq(-abs(max_bk),0,by=abs(max_bk)/50),seq(0,abs(max_bk),by=abs(max_bk)/50))) # 为确保0值为白色
  p <- NA
  
  
  # 处理只有一个代谢物的情况
  if(ncol(da_plot2) == 1 | nrow(da_plot2) == 1){
    p <- pheatmap::pheatmap(as.matrix(da_plot2), scale = "none",  legend = T, angle_col = "45", 
                            cluster_cols = F, cluster_rows = F,
                            width = 20, height = 10,cellwidth = 20, cellheight = 30,breaks = bk,
                            color = c(colorRampPalette(colors = c("#30A9DE","white"))(length(bk)/2),
                                      colorRampPalette(colors = c("white","#E53A40"))(length(bk)/2)))
    mat <- as.matrix(da_plot2)
  }else{
    p <- pheatmap::pheatmap(as.matrix(da_plot2), scale = "none",  legend = T, angle_col = "45", 
                            width = 20, height = 10,cellwidth = 20, cellheight = 30,breaks = bk,
                            color = c(colorRampPalette(colors = c("#30A9DE","white"))(length(bk)/2),
                                      colorRampPalette(colors = c("white","#E53A40"))(length(bk)/2)))
    row_order <- p$tree_row$order
    ordered_rownames <- rownames(da_plot2)[row_order]
    col_order <- p$tree_col$order
    ordered_colnames <- colnames(da_plot2)[col_order]
    mat <- as.matrix(da_plot2[ordered_rownames[length(ordered_rownames):1],ordered_colnames])
    
  }
  
  
  if(database %in% c("Class","Pathway","Disease","Drug") ){
    fn=paste("HUBMet",database, sep="_")
  }else{
    fn=database
  }
  
  
  
  pdf(paste("./Output/M2_MSEA/",fn,"_heatmap.pdf",sep = ""),width = 11,height = 8)
  print(p) 
  dev.off()
  
   
  
}







## --------------------Running score plot for MSEA----------------

 
msea_ER_forplot <- function(term_met, met_id, met_id_fc, forplot = F){
  # 该函数是对单个代谢物集合进行单个metabolite set的enrichment分析
  # 最终输出结果为最高的ES, pvalue，NES
  # 以及每个数据集的N次permutation的NES结果，用于计算FDR
  
  
  term_met <- strsplit(term_met, split = ";")[[1]]
  overlap <- intersect(term_met, met_id)
  
  # 设定permutation的次数， 1000次,代码测试时可设置较小数值
  #permutation_times <- 1000
  
  if(length(overlap) == 0){ # 如果input和数据集没有overlap，则输出NA
    if(forplot == T){
      return(NA)
    }else{
      return(c(0,NA,NA, NA, NA, rep(NA, times = permutation_times)))
    }
    
  }else{
    # 排序，check overlap
    da <- data.frame(met_id, met_id_fc)
    da$met_id_fc <- as.numeric(da$met_id_fc)
    
    da <- solve_redundant(da$met_id, da$met_id_fc, choice="mean")
    colnames(da) <- c("met_id","met_id_fc")
    da <- da[order(da$met_id_fc, decreasing = T),] # 按照FC排序
    da$rank <- 1:nrow(da)
    rownames(da) <- da$met_id
    
    da[which(da$met_id %in% overlap),"hit"] <- "yes"
    da[which(!(da$met_id %in% overlap)),"hit"] <- "no"
    
    # 计算real enrichment score
    sum_hit <- sum(abs(da[which(da$hit == "yes"),"met_id_fc"]))  # sum of abs
    da[which(da$hit == "yes"),"score"] <- abs(da[which(da$hit == "yes"),"met_id_fc"])/sum_hit # fraction of abs
    da[which(da$hit == "no"),"score"] <- -1/(length(met_id)-length(overlap))
    
    
    for(i in 1:nrow(da)){
      if(i == 1){
        da[i,"ES"] <- da[i, "score"]
      }else{
        da[i,"ES"] <- da[i, "score"] + da[i-1,"ES"]
      }
    }
    da$ES_abs <- abs(da$ES)
    da[which(da$ES_abs == max(da$ES_abs)),"label"] = "max"
    return(da)
  }
}


msea_plot_ES_bar <- function(database, pathway_ID, hmdb_list, hmdb_list_value, res_msea){
  
  
  library(tidyverse)
  choice = "mean"
  # 处理res为空data frame的情况
  if(nrow(res_msea) == 0){
    p <- ggplot()+theme(panel.background = element_rect(fill='transparent'))
    return(p)
  }
  
  res <- res_msea
  rm(res_msea)
  
  
  # 获取pathway信息
  if(database == "smpdb"){
    term_met <- smpdb_term_anno_new_MSEA[which(smpdb_term_anno_new_MSEA$SMP.ID == pathway_ID),"HMDB.ID"]
    path_name <- smpdb_term_anno_new_MSEA[which(smpdb_term_anno_new_MSEA$SMP.ID == pathway_ID),"Pathway.Name"]
    nsp <- res[which(res$SMP.ID == pathway_ID), c("ES","NES","FDR","pvalue_nominal")]
    n_check <- nrow(res[which(res$SMP.ID == pathway_ID), ])
  }else if(database == "kegg"){
    term_met <- kegg_term_anno_new[which(kegg_term_anno_new$KEGG_pathway_ID == pathway_ID),"HMDB.ID"]
    path_name <- kegg_term_anno_new[which(kegg_term_anno_new$KEGG_pathway_ID == pathway_ID),"KEGG_pathway_name"]
    nsp <- res[which(res$KEGG_pathway_ID == pathway_ID), c("ES","NES","FDR","pvalue_nominal")]
    n_check <- nrow(res[which(res$KEGG_pathway_ID == pathway_ID), ])
  }else if(database == "humanGEM"){
    term_met <- humanGEM_term_new[which(humanGEM_term_new$SUBSYSTEM == pathway_ID),"HMDB.ID"]
    path_name <- pathway_ID
    nsp <- res[which(res$SUBSYSTEM == pathway_ID), c("ES","NES","FDR","pvalue_nominal")]
    n_check <- nrow(res[which(res$SUBSYSTEM == pathway_ID), ])
    
  }else if(database == "reactome"){
    term_met <- reactome_term_anno_new[which(reactome_term_anno_new$reactome_id == pathway_ID),"HMDB.ID"]
    path_name <- reactome_term_anno_new[which(reactome_term_anno_new$reactome_id == pathway_ID),"reactome_name"]
    nsp <- res[which(res$reactome_id == pathway_ID), c("ES","NES","FDR","pvalue_nominal")]
    n_check <- nrow(res[which(res$reactome_id == pathway_ID), ])
  }else if(database == "disease"){
    term_met <- disease_term_new[which(disease_term_new$disease == pathway_ID),"HMDB.ID"]
    path_name <- pathway_ID
    nsp <- res[which(res$disease == pathway_ID), c("ES","NES","FDR","pvalue_nominal")]
    n_check <- nrow(res[which(res$disease == pathway_ID), ])
  }else if(database %in% c("Class","Pathway","Disease","Drug")){ # 0322, HBM适应化处理
    term_met <- HUBMet_term[which(HUBMet_term$Term == pathway_ID),"HBM_ID"]
    path_name <- pathway_ID
    nsp <- res[which(res$Term == pathway_ID), c("ES","NES","FDR","pvalue_nominal")]
    n_check <- nrow(res[which(res$Term == pathway_ID), ])
  }
  
  
  # 计算画图所需running score
  
  if(n_check==0){
    print("No hit metabolite in the metabolite sets! or No enough metabolite in the data set!")
  }else{
    temp <- solve_redundant(hmdb_list , hmdb_list_value, choice)
    es_path <- msea_ER_forplot(term_met, temp$hmdb_list , temp$hmdb_list_value, forplot = T) # 1103
    es <- es_path[which(es_path$label == "max"), "ES"]
    
    #y_label = max(es_path$ES) - 0.1*(abs(min(es_path$ES))+abs(max(es_path$ES))) #  remove
    # 修改label_ES
    label_ES <- paste( c(paste("ES = ", round(nsp[1],3), sep = ""), 
                         paste("NES = ",round(nsp[2],3), sep = ""), 
                         paste("P = ",round(nsp[4],3), sep = ""), 
                         paste("FDR = ", round(nsp[3],3), sep = "")), collapse  = "; ") 
    
    # pa有较多修改
    pa <- ggplot(es_path)+
      geom_hline(aes(yintercept=0),linetype="dashed",color="black", linewidth=0.5)+
      geom_line( mapping = aes(x=rank, y = ES), color = "grey20")+
      geom_bar(subset(es_path, label == "max" ) , mapping = aes(x=rank, y = ES), 
               
               color = "red",fill = "red",stat = "identity", width = 0.25)+ #0813 添加fill = "red",
      #scale_x_continuous(limits = c(1, max(es_path$rank)) )+
      scale_x_continuous(limits = c(0, max(es_path$rank)+1) )+ # 0813 修改
      labs(x=element_blank(), y="Running Enrichment Score"   )+   
      ggtitle(paste( unique(c(path_name, pathway_ID)), collapse = "\n"), label_ES) +  
      
      theme_bw() +
      scale_y_continuous(expand = c(0,0))+
      theme(panel.grid=element_blank(),
            #axis.line = element_line(linewidth = 0.5),
            plot.title = element_text(hjust = 0.5,size=14, colour = "black", face = "bold"),
            plot.subtitle = element_text(hjust = 0.5,size=13, colour = "black", face="italic"),  
            panel.background = element_blank(),
            panel.border = element_rect(linewidth = 0.5, color="black"),
            #strip.text = element_text(face = "bold", size = 13),
            legend.position = "none",
            axis.ticks.x = element_blank(),
            axis.title.x =element_blank(), 
            axis.text.x =element_blank(),
            plot.margin = unit(c(0.2,0.2,0,0.2), "cm"),
            axis.title.y =element_text(size=13, colour = "black"),
            axis.text.y = element_text(size=13, colour = "black"))
    
    # 增加pm
    pm <- ggplot(es_path)+
      #scale_x_continuous(limits = c(1, max(es_path$rank)) )+
      scale_x_continuous(limits = c(0, max(es_path$rank)+1) )+ # 0813 修改
      geom_vline(xintercept = es_path[which(es_path$hit == "yes"),]$rank,linetype="solid",color="black", linewidth=0.5)+
      theme_bw()+
      theme(axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.title.x = element_blank(), # 0120添加
            panel.border = element_rect(linewidth = 0.5, color="black"),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank())
    
    
    pb <- ggplot(es_path)+
      
      geom_bar(subset(es_path, met_id_fc>=0) ,mapping = aes(x=rank, y = met_id_fc),
               color = "#E53A40", fill= "#E53A40", stat = "identity", width = 0.1)+  # 添加 fill= "#E53A40",
      geom_bar(subset(es_path,met_id_fc<0) ,mapping = aes(x=rank, y = met_id_fc),
               color = "#30A9DE", fill= "#30A9DE", stat = "identity", width = 0.1)+ # 添加 fill= "#30A9DE",
      geom_hline(aes(yintercept=0),linetype=5,col="black", linewidth=0.5)+
      #scale_x_continuous(limits = c(1, max(es_path$rank)) )+
      scale_x_continuous(limits = c(0, max(es_path$rank)+1) )+ # 0813 修改
      labs(x="Metabolite rank", y="Ranking metric", title = element_blank())+
      theme_bw() +
      theme(panel.grid=element_blank(),
            panel.border = element_rect(linewidth = 0.5, color="black"),
            #plot.title = element_blank(),
            panel.background = element_blank(),
            #strip.text = element_text(face = "bold", size = 12),
            legend.position = "none",
            plot.margin = unit(c(0.0,0.2,0.2,0.2), "cm"),
            axis.title  =element_text(size=13, colour = "black"), 
            #axis.text.x =element_text(size=10, colour = "black"),
            #axis.title.y =element_text(size=12, colour = "black"),
            axis.text = element_text(size=13, colour = "black"))  
    
    suppressWarnings({p <- cowplot::plot_grid(pa,pm,NULL,pb,ncol = 1, align = "v", rel_heights = c(1,0.2,-0.0225,0.3))}) 
    # add suppressWarnings for 'remove NS'
    
    # 控制某些时候pathway_ID中存在"/"的情况，会被误识别为路径符
    pathway_ID <- gsub("/","_",pathway_ID)
    
    if(database %in% c("Class","Pathway","Disease","Drug")){
      database <- paste("HUBMet_",database,sep="")
    }
    
    
    ggsave(filename = paste("./Output/M2_MSEA/",database,"_",pathway_ID,"_RES.pdf",sep = ""),p,width = 7,height = 7) 
    #return(p)
  }
  
}



# ------------------ Pathview in ORA and MSEA---------------------


#  此函数只针对kegg smpdb reactome 这3个数据库进行设置
# 在结果显示页面提供单独的下载链接

pathview_kegg_smpdb_reactome <- function(database, pathway_ID, hmdb_list, hmdb_list_value = NA){
  library(pathview)
  library(SBGNview)
  library(SBGNview.data)
  library(tidyverse)
  data(sbgn.xmls) 
  
 
  dir.create("./Output/pathview", recursive = TRUE, showWarnings = FALSE)
  if(sum(!is.na(hmdb_list_value))==length(hmdb_list)){
    da <- data.frame("HMDB_ID_unique"=hmdb_list,
                     "cpd_da" = hmdb_list_value)
    data_flag = T
  }else{
    da <- data.frame("HMDB_ID_unique"=hmdb_list,
                     "cpd_da" = rep(0, times = length(hmdb_list)))
    data_flag = F
  }
  
  da <- left_join(da, HMDB_ID_unique_name[,c("HMDB_ID_unique","KEGG_Compound_ID")],by = join_by(HMDB_ID_unique))
  da <- da[which(!is.na(da$KEGG_Compound_ID)),]
  
  # 处理KEGG ID冗余问题, 画图默认取mean
  if(length(da$KEGG_Compound_ID) > length(unique(da$KEGG_Compound_ID))){
    a <- data.frame(table(da$KEGG_Compound_ID))
    a$Var1 <- as.character(a$Var1)
    a <- a[which(a$Freq>1),]
    for(i in a$Var1){
      da[which(da$KEGG_Compound_ID == i),"cpd_da" ] <- mean(da[which(da$KEGG_Compound_ID == i),"cpd_da" ])
    }
    da <- unique(da[,c(2,3)])
  }
  
  if(database == "kegg"){
    met_list_keggID <- da$KEGG_Compound_ID
    cpd_da <- da$cpd_da
    names(cpd_da) <- met_list_keggID
    p_id <- gsub(x = pathway_ID, replacement = "" , pattern = "map")
    if(data_flag == F){
      pathview(cpd.data = cpd_da, cpd.idtype = "kegg", kegg.native =T,map.symbol = TRUE, map.cpdname = TRUE,
               pathway.id = p_id, res=300,
               kegg.legend = ("node"),
               kegg.dir ="./Output/pathview/", same.layer = T,
               mid = list(cpd="yellow"),low=list(cpd="white"),high=list(cpd="white"),
               species = "hsa", out.suffix ="need")
    }else{
      pathview(cpd.data = cpd_da, cpd.idtype = "kegg", kegg.native =T,map.symbol = TRUE, map.cpdname = TRUE,
               pathway.id = p_id, res=300,
               kegg.legend = ("node"),
               kegg.dir ="./Output/pathview/", same.layer = T,
               mid = list(cpd="grey40"),low=list(cpd="#30A9DE"),high=list(cpd="#E53A40"),
               species = "hsa", out.suffix ="need")
    }
    
    file.remove(paste(c("./Output/pathview/","hsa",p_id,".png"), collapse = ""))
    file.remove(paste(c("./Output/pathview/","hsa",p_id,".xml"), collapse = ""))
    file.rename(paste(c("./","hsa",p_id,".need.png"), collapse = ""), 
                paste(c("./Output/pathview/","hsa",p_id,".need.png"), collapse = ""))
    filename <- paste(c("./Output/pathview/","hsa",p_id,".need.png"), collapse = "")
    # 将pathview直接展现出来
    img <- png::readPNG(filename, native =T, info = T)
    grid::grid.newpage()
    grid::grid.raster(img)
    
  }else{
    if(database == "smpdb"){
      pathway_ID <- gsub(pattern = "SMP00", replacement = "SMP", pathway_ID)
    }
    data(pathways.info)
    if(pathway_ID %in% pathways.info$pathway.id){
      met_list_keggID <- da$KEGG_Compound_ID
      cpd_da <- da$cpd_da
      
      cpd_da <- data.frame("value" = cpd_da)
      rownames(cpd_da) <- met_list_keggID
      if( data_flag == F ){
        SBGNview.obj <- SBGNview(cpd.data =cpd_da,
                                 cpd.id.type = "kegg",
                                 input.sbgn = pathway_ID ,
                                 col.cpd.low = "white",
                                 col.cpd.high = "white",
                                 col.cpd.mid = "yellow",
                                 output.file = "./Output/pathview/need", 
                                 output.formats =  c("png")) 
      }else{
        SBGNview.obj <- SBGNview(cpd.data =cpd_da,
                                 cpd.id.type = "kegg",
                                 input.sbgn = pathway_ID ,
                                 col.cpd.low = "#30A9DE",
                                 col.cpd.high = "#E53A40",
                                 col.cpd.mid = "grey40",
                                 output.file = "./Output/pathview/need", 
                                 output.formats =  c("png")) 
      }
      
      
      print(SBGNview.obj) # 输出png图片
      
      
      if(database == "smpdb"){
        file.remove(paste(c("./http___identifiers.org_","smpdb_",pathway_ID,".sbgn"), collapse = ""))
      }else{
        file.remove(paste(c("./http___identifiers.org_","reactome_",pathway_ID,".sbgn"), collapse = ""))
      }
      
      file.remove(paste(c("./Output/pathview/","need_",pathway_ID,".svg"), collapse = ""))
      filename <- paste(c("./Output/pathview/","need_",pathway_ID,".png"), collapse = "")
      img <- png::readPNG(filename, native =T)
      grid::grid.newpage()
      grid::grid.raster(img)
      
    }else{
      return("PathView not avaliable!")
    }
  }
  
  
}


