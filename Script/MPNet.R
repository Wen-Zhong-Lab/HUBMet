# Functional module 4: MPNet

# -------------------Load data and library -------------

library(plyr)
library(tidyverse)
library(igraph)
library(qgraph)
 
dir.create("./Output/M4_MPNet", recursive = TRUE, showWarnings = FALSE)

HUBMet_annotation <- readRDS("./Data/HUBMet_annotation.RDS") 
pro_met_ass_anno <- readRDS("./Data/pro_met_ass_anno.RDS")
pro_met_ass <- pro_met_ass_anno[,c("PID","HBM_ID","Index","Evidence")]
proID_ref_annotation <- readRDS("./Data/proID_ref_annotation.RDS") 


# ------------- Main function ------------------------------

ProMetNetwork <- function(proteins = NA, metabolites = NA, 
                          protein_type = "Uniprot", job_id ="jobid",
                          timen = 0){
  library(qgraph)
  library(igraph) 
  library(ggsci)  
  evidence = "Validated"  
  # limit the range of expanding
  if(timen > 5){
    timen <- 5
  }
  if(timen < 0){
    timen <- 0
  }
  
  
  # needed data: pro_met_ass_anno, proID_ref_annotation, HUBMet_annotation
  # Output: one link table, one node table, 3 bar plots of community characteristics, one network in the folder 'Output/M4_MPNet'
  
  flag_noPro = NULL
  flag_noMet = NULL
  
  
  if( "all" %in% evidence  ){
    pro_met_ass_anno_used <- pro_met_ass_anno
  }else{
    pro_met_ass_anno_used <- pro_met_ass_anno %>%
      filter(Evidence %in% setdiff(evidence,"all"))
  }
  
  
  # data preparation
  
  # check proteins
  temp <- proID_ref_annotation
  colnames(temp)[colnames(temp) == protein_type] <- "Input"
  temp <- na.omit(unique(temp[,c("PID","Input")]))
  pid_input <- na.omit(unique(temp[which(temp$Input %in% na.omit(unique(proteins))),]))
  pid <- na.omit(unique(temp[which(temp$Input %in% na.omit(unique(proteins))),]$PID))
  
  pid <- na.omit(intersect(pid, pro_met_ass_anno_used$PID))
  pid_input <- pid_input[which(pid_input$PID %in% pid),]
  
  if(length(pid) == 0){
    flag_noPro = TRUE
  }else{
    flag_noPro = FALSE
  }
  
  
  # check metabolite
  metabolites <- unique(intersect( na.omit(metabolites[!(metabolites %in% c("","NA"))]), pro_met_ass_anno_used$HBM_ID))
  if(length(metabolites) == 0){
    flag_noMet = TRUE
  }else{
    flag_noMet = FALSE
  }
  
  # section to deal with different input situation
  if(flag_noPro == TRUE & flag_noMet == TRUE){
    print("There is no protein and metabolite mapped in the Protein-Metabolite Network Database.")
    return("Analysis terminated.")
    
  }else if(flag_noPro == TRUE & flag_noMet == FALSE){
    print("No inputted protein was mapped to the Protein-Metabolite Network Database Or no protein was inputted.")
    print("The protein(s) connected to the inputted metabolite(s) was used for analysis.")
    link <- pro_met_ass_anno_used[which(pro_met_ass_anno_used$HBM_ID %in% metabolites),]
    
    if(nrow(link) == 0){
      print("But no link between protein and metabolite was found!")
      return("Analysis terminated.")
    }
    
    pid <- na.omit(unique(link$PID))
    
  }else if(flag_noPro == FALSE & flag_noMet == TRUE){
    print("No inputted metabolite was mapped to the Protein-Metabolite Network Database Or no metabolite was inputted.")
    print("The metabolite(s) connected to the inputted protein(s) was used for analysis.")
    link <- pro_met_ass_anno_used[which(pro_met_ass_anno_used$PID %in% pid),]
    
    if(nrow(link) == 0){
      print("But no link between protein and metabolite was found!")
      return("Analysis terminated.")
    }
    
    metabolites <- na.omit(unique(link$HBM_ID))
    
  }else{  
    link <- pro_met_ass_anno_used[which(pro_met_ass_anno_used$HBM_ID %in% metabolites & pro_met_ass_anno_used$PID %in% pid),]
    if(nrow(link) == 0){
      
      print(paste(c(length(unique(link$PID))," protein(s) and ",
                    length(unique(link$HBM_ID)),
                    " metabolite(s) were matched in the network database. But no link was found between them."), collapse =""))
      
      return("Analysis terminated.")
    }
  }
  
  
  
  
  # when there are links detected, keep on the following code
  node <- data.frame("name" = append( unique(link$PID), unique(link$HBM_ID)),
                     "type" = append( rep("protein", times = length(unique(link$PID))) , 
                                      rep("metabolite", times = length(unique(unique(link$HBM_ID)))) ))
  
  # Print the result of database mapping
  
  print(paste(c(length(unique(link$PID))," protein(s) and ",
                length(unique(link$HBM_ID)),
                " metabolite(s) were matched in the network database. ",
                nrow(link), " link(s) were found between them."), collapse =""))
  
  # network construction
  network <- igraph::graph_from_data_frame(d=na.omit(link[,c("PID","HBM_ID")]), 
                                           directed=F,  vertices = na.omit(node[,c("name","type")])) 
  
  # parameter calculation
  rownames(node) <- node$name
  temp <- as.data.frame( igraph::degree(network))
  colnames(temp)[1] <- "degree"
  node <- cbind(node, temp)
  temp <- as.data.frame( igraph::closeness(network))
  colnames(temp)[1] <- "closeness"
  node <- cbind(node, temp)
  temp <- as.data.frame( igraph::betweenness(network))
  colnames(temp)[1] <- "betweenness"
  node <- cbind(node, temp)
  
  # community detection
  set.seed(1) 
  temp <- as.data.frame( igraph::membership(igraph::cluster_fast_greedy(network))) 
  colnames(temp)[1] <- "community"
  node <- cbind(node, temp)
  
  # add necessary annotation
  node <- dplyr::left_join(node, unique(na.omit(HUBMet_annotation[,c("HBM_ID","Name","Category")])),
                           by=c("name"="HBM_ID")) 
  
  node[which(node$name == "HBM00001818" ),"Name"] <- "NAD+"
  node[which(node$name == "HBM00002633" ),"Name"] <- "ATP"
  node[which(node$name == "HBM00001025" ),"Name"] <- "ADP"
  node[which(node$name == "HBM00001336" ),"Name"] <- "AMP"
  
  
  node <- dplyr::left_join(node, unique(na.omit(proID_ref_annotation[,c("PID","Label","Function")])), by=c("name"="PID"))
  node[which(is.na(node$Name)),"Name"] <- node[which(is.na(node$Name)),"Label"]
  node[which(is.na(node$Category)),"Category"] <- node[which(is.na(node$Category)),"Function"]
  node <- unique(node[,1:8])
  
  colnames(node)[c(1,2,7)] <- c("id","type","name")
  node <- node[,c(1,2,7,3:6,8)]
  
  # prepare results for export
  node <- node[order(-node$degree,-node$closeness),]
  node_export <- left_join(node, pid_input,by=c("id"="PID"))
  colnames(node_export)[which(colnames(node_export) == "Input")] <- "Protein.Input"   
  node_export <- node_export[,c(1:3,9,8,4:7)]
  node_export[is.na(node_export)] <- ""
  colnames(node_export) <- c("ID","Type","Name","Protein.Input","Class.Function",   
                             "Degree","Closeness","Betweenness","Community")
  
  link_export <- link[,c("PID","HBM_ID","Evidence")]
  colnames(link_export)[c(1,2)] <- c("Protein.ID","HUBMet.ID") # 1029 add and change
  link_export[is.na(link_export)] <- ""
  write.table(link_export, paste("./Output/M4_MPNet/link_",job_id,".txt",sep="") ,
              sep = "\t", quote = F, row.names = F)
  write.table(node_export, paste("./Output/M4_MPNet/node_",job_id,".txt",sep="") ,
              sep = "\t", quote = F, row.names = F)
  
  
  theme_basic = theme_classic() +
    theme(panel.grid=element_blank(),panel.background = element_blank(),
          legend.position = "top",
          axis.line = element_line(color="black", linewidth = 0.5),
          axis.title  =element_text(size=12, colour = "black"),  
          axis.text  =element_text(size=12, colour = "black"),  
          legend.text=element_text(size=12, colour = "black"), 
          legend.title=element_text(size=12, colour = "black"))
  
  # figure 1 bar plot of node number in community

   
  p1 <- na.omit(unique(node[,c("id","type","community")])) %>% 
    mutate(gp = paste0(type,"_", community ))
  
  p1 <- as.data.frame(table(p1$gp))
  p1$Var1 <- as.character(p1$Var1)
  p1$type="protein"
  p1[which(grepl("metabolite",p1$Var1)),"type"] <- "metabolite"
  p1$community <- as.integer(unlist(lapply(p1$Var1,function(x){strsplit(x,"_")[[1]][2]})))
  p1$n <- p1$Freq
  
  
  com_include <- unique(p1[which(p1$n > 5),]$community)
  
  if(length(com_include) == 0){
    com_include <- unique(p1$community)
  }
  
  p1 <- p1 %>% filter(community %in% com_include) 
  p1x <- p1 %>% 
    mutate(type=factor(type,levels=c("metabolite","protein")[2:1])) %>% 
    ggplot(aes(as.factor(community),n))+
    geom_bar(aes(fill=type), stat="identity", position="stack")+
    scale_y_continuous(expand = c(0,0), limits = c(0,max(p1$n)*1.05 ))+
    scale_fill_manual(values = c("metabolite"="#C6E3A1", "protein" ="#FED7D9"))+
    labs(x="Community", y="Number of metabolites and proteins",fill="Type")+
    theme_basic
  
  
  ggsave(paste("./Output/M4_MPNet/network_",job_id,"_figure1.pdf",sep = ""),p1x, width = 4.125, height = 4)
  
  
   
  
  # figure 2
  temp <- unique(node[which(node$type == "metabolite"),c("id","Category")])
  count_class <- as.data.frame(table(temp$Category))
  
  node_community <- node[which(node$type == "metabolite"),]
  ncom_met <- unique( node_community$community)   
  for(i in ncom_met){   
    temp <- node_community[which(node_community$community == i ),]
    a <- as.data.frame(table(temp$Category))
    a <- left_join(a, count_class,by="Var1") 
    a$perc <- 100*(a$Freq.x/a$Freq.y)
    a$com <- i
    a$n_com <- length(unique(temp$id))
    if(i == ncom_met[1]){   
      need <- a
    }else{
      need <- rbind(need,a)
    }
  }
  
  class9_level <- c("Lipid","Xenobiotics","Amino Acid","Peptide","Carbohydrate",
                    "Nucleotide","Cofactors and Vitamins","Energy","Others") 
  
  color_class9 <- c('#ffa510','#45aee8','#EB99C2','#f6c9aa','#f57070','#add547','#41b7ac','#ad85d5','#C3DEE0')
  names(color_class9) <- class9_level
  color_class9_df <- as.data.frame(color_class9) %>% 
    rownames_to_column(var="Category")
  
  
  need$Var1 <- factor(need$Var1, levels = class9_level[9:1])
  need$label <- paste(round(need$perc,0),"%",sep = "")
  
  need$com_text <- paste( "C.", need$com, "\nn=",need$n_com,sep = "")
  need$text <- paste(need$Var1, "\nn=",need$Freq.y,sep="")
  
  xc <- unique(need[,c("com","com_text")]) %>%
    arrange(com)
  need$com_text <- factor(need$com_text, levels = xc$com_text)
  
  xn <- unique(need[,c("Var1","text")])
  
  yorder_met <- unique(need[,c("Var1","Freq.y")]) %>% 
    mutate(Var1 = as.character(Var1)) %>% 
    filter(Var1 %in% as.character(unique(need[which(need$com %in% com_include),]$Var1))) %>% 
    arrange(-Freq.y)
  yorder_met <- yorder_met$Var1
  
  need <- need[which(need$com %in% com_include),]
  need$Var1 <- factor(as.character(need$Var1), levels = yorder_met[length(yorder_met):1])
  
  
  p2 <- ggplot(need, aes(x=Var1,y=perc,fill=Var1))+
    geom_bar(stat = "identity" )+
    geom_text(mapping = aes(x=Var1,y=0,label = label),hjust=0,vjust=0.5, size=3)+
    labs(y="Percentage of metabolites in each community",x="Category of metabolite")+
    scale_x_discrete(breaks = xn$Var1, label=xn$text)+
    coord_flip() + 
    scale_fill_manual(values = color_class9)+
    theme_basic+
    theme(legend.position = "none",
          axis.text.y = element_text(size = 10),  
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())+ 
    facet_wrap(~com_text,nrow = 1)
  ggsave(paste("./Output/M4_MPNet/network_",job_id,"_figure2.pdf",sep = ""),p2, width = 8.25, height = 4)
  
  
   
  
  
  # figure 3
  da_com <- node[which(node$type == "protein"),]
  da_anno <- proID_ref_annotation %>% 
    filter(PID %in% da_com$id) %>%
    select(PID, HPA_protein_class) %>%
    unique()
  da_anno[is.na(da_anno$HPA_protein_class),2] <- "Unknown"
  
  com_n <- as.data.frame(table(da_com$community)) %>%
    mutate(Var1=as.numeric(Var1))
  
  need <- do.call(rbind,apply(unique(da_anno), 1, function(x){return(data.frame("PID"=x[1],
                                                                                "class"=strsplit(x[2],";")[[1]],
                                                                                row.names = NULL))})) %>%
    unique() %>%
    left_join(da_com[,c("id","community")], by=c("PID"="id"))
  
  hpa_class <- need[,1:2]
  
  hpa_class <- as.data.frame(table(hpa_class$class)) %>% arrange(-Freq)
  
  hpa_class$need <- "no"
  hpa_class[1:10,"need"] <- "yes"
  hpa_class <- na.omit(hpa_class)
  hpa_class$class <- factor(hpa_class$Var1, levels = hpa_class$Var1)
  hpa_class$n <- hpa_class$Freq
  
  
  p1 <- need %>% 
    mutate(gp = paste0(class,"_", community ))
  
  p1 <- as.data.frame(table(p1$gp))
  p1$Var1 <- as.character(p1$Var1)
  p1$class <- as.character(unlist(lapply(p1$Var1,function(x){strsplit(x,"_")[[1]][1]})))
  p1$community <- as.integer(unlist(lapply(p1$Var1,function(x){strsplit(x,"_")[[1]][2]})))
  p1$n <- p1$Freq
  
  
  
  daforplot <- p1[,c("class","community","n")] %>%
     left_join(hpa_class[,c("class","n")], by="class") %>%
    mutate(perc = paste(round(100*n.x/n.y,0), "%",sep=""),
           perc_n = 100*n.x/n.y  ) %>%
    left_join(com_n, by=c("community"="Var1")) %>%  
    mutate(com_n = paste( "C.", community,"\nn=", Freq,sep=""))
  
  dac <- unique(daforplot[,c("community","com_n")]) %>%
    arrange(community)  
  daforplot$com_n <- factor(daforplot$com_n, levels = dac$com_n)
  
  class_need <- c()
  for(i in unique(daforplot$community)){   
    temp <- daforplot[which(daforplot$community == i),] %>%
      arrange(-n.x, -perc_n)
    class_need <- c(class_need, temp$class[1:3])
    
  }
  class_need <- setdiff(class_need, c("NA",NA,"Unknown"))
  
  
  xn <- unique(daforplot[,c("class", "n.y")]) %>%
    arrange(-n.y) %>%
    mutate(text = paste(  class , "\nn=",n.y, sep=""))
  
  daforplot$class<- factor(daforplot$class, levels=xn$class[nrow(xn):1])
  
  
  daforplot <- daforplot[which(daforplot$class %in% class_need & 
                                 daforplot$community %in% com_include),]
  
  orderx <- daforplot[,c("class","n.y")] %>% unique() %>% arrange(-n.y)
  
  orderx <- as.character(orderx$class)
  
  color_proteinclass <- pal_d3("category20")(length(orderx) )
  names(color_proteinclass) <- orderx
  color_proteinclass_df <- as.data.frame(color_proteinclass) %>% 
    rownames_to_column(var="class")
  
  
  p3 <- ggplot(daforplot, aes(x=class,y=perc_n,fill=class))+
    geom_bar(stat = "identity" ,alpha=0.75 )+
    geom_text(mapping = aes(x=class,y=0,label = perc),hjust=0,vjust=0.5, size=3)+
    labs(y="Percentage of proteins in each community",x="Protein Class\nHuman Protein Atlas")+
    scale_x_discrete(breaks = xn$class, label=xn$text)+
    coord_flip() + 
    scale_fill_manual(values = color_proteinclass)+
    theme_basic+
    theme(legend.position = "none",
          axis.text.y = element_text(size = 10 ),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())+ 
    facet_wrap(~com_n,nrow = 1)
  ggsave(paste("./Output/M4_MPNet/network_",job_id,"_figure3.pdf",sep = ""),p3, width = 8.25, height = 4)
  
 
  
  # chose part of the network to visualize, depends on the number of link
  link_limits <- 5*(10+timen*2)
  node_limiks <- 10+timen*2
  if( nrow(link) >= link_limits ){   
    node <- node[which(!(node$degree == 1 & node$closeness == 1)),]
    temp1 <- node[which(node$type == "protein"),]
    temp1 <- temp1[order(-temp1$degree,-temp1$closeness),]
    temp1 <- na.omit(temp1[1:node_limiks,]) 
    
    temp2 <- node[which(node$type == "metabolite"),]
    temp2 <- temp2[order( -temp2$degree,-temp2$closeness),]
    temp2 <- na.omit(temp2[1:node_limiks,])  
    need_index <- c("PID","HBM_ID","edge_weight")
    linkPlot <- link[which(link$HBM_ID %in% temp2$id | link$PID %in% temp1$id),]
    
    
    if(nrow(linkPlot) >= link_limits){
      # if it is still larger than 50, control the number of links to be around 50
      linkPlot_part1 <- link[which(link$HBM_ID %in% temp2$id & link$PID %in% temp1$id), ]
      if(nrow(linkPlot_part1) <= link_limits){
        nlink <- link_limits-nrow(linkPlot_part1)
        linkPlot_part2 <- link[which(link$HBM_ID %in% temp2$id), ][1:ceiling(nlink/2),]
        linkPlot_part3 <- link[which(link$PID %in% temp1$id), ][1:ceiling(nlink/2),]
        linkPlot <- unique(rbind(rbind(linkPlot_part1, linkPlot_part2), linkPlot_part3))
        linkPlot <- na.omit(linkPlot)  
        rm(linkPlot_part1)
        rm(linkPlot_part2)
        rm(linkPlot_part3)
      }else{
        linkPlot <- linkPlot_part1
      }
    }
    
    # construct the network again based on filtered links
    nodePlot <- node[which(node$id %in% append(linkPlot$PID, linkPlot$HBM_ID)),]
    
    network <- igraph::graph_from_data_frame(d=na.omit(linkPlot[,c("PID","HBM_ID","edge_weight")]), 
                                             directed=F,  vertices =na.omit(nodePlot[,c("id","type")]))
    
    plot_size <- max(8, 8*nrow(linkPlot)/50)  
    # control the size of node: 5-30
    daforsize <- data.frame("degree"=unique(nodePlot$degree),"size"=NA)
    if(max(daforsize$degree) <= 30 & min(daforsize$degree) >= 5){
      V(network)$size <- nodePlot$degree
    }else{
      step = (30-5)/(nrow(daforsize)-1)
      daforsize <- daforsize[order(-daforsize$degree),]
      daforsize$size <- seq(30,5,-step)[1:nrow(daforsize)]
      nodePlot <- dplyr::left_join(nodePlot, daforsize,by="degree")
      V(network)$size <- nodePlot$size
      nodePlot <- nodePlot[,which(colnames(nodePlot) != "size")]
    }
    
    V(network)$name <- nodePlot$name
    V(network)$color <- ifelse(V(network)$type == "metabolite", "#C6E3A1", "#FED7D9")
    
    # collect community information for visualization
    com <- list()
    for(i in unique(nodePlot$community)){
      com[[paste("com",i,sep = "")]] <- nodePlot[which(nodePlot$community == i),]$name
    }
    
    
    e <- igraph::as_edgelist(network,names=FALSE)
    l <- qgraph::qgraph.layout.fruchtermanreingold(e, vcount = vcount(network),
                                                   area = 6*(vcount(network)^2), repulse.rad = vcount(network)^3)
    
    
    # export network plot
    pdf(paste("./Output/M4_MPNet/network_",job_id,"_time_",timen,".pdf",sep = ""), width = plot_size, height = plot_size)
    set.seed(1)
    plot.igraph(network, 
                layout = l,
                mark.groups = com,
                mark.col = rgb(0.8,0.8,0.8,0.25), 
                mark.border = "grey85", 
                vertex.lable =V(network)$name,
                vertex.color = V(network)$color,
                vertex.size = V(network)$size,
                vertex.label.color = "black",
                vertex.label.family = "sans",  
                edge.color = "grey45", 
                edge.width=linkPlot$edge_weight*1.5,
                vertex.frame.color ="white",
                margin = rep(0, times = 4) )
    dev.off()
  }else{
    # if the number of links is less than 50, visualize all the nodes and links
    linkPlot <- link
    nodePlot <- node
    
    network <- igraph::graph_from_data_frame(d=na.omit(linkPlot[,c("PID","HBM_ID","edge_weight")]), 
                                             directed=F,  vertices =na.omit(nodePlot[,c("id","type")]))
    plot_size <- max(8, 8*nrow(linkPlot)/50)  
    if(nrow(linkPlot) == 1){
      V(network)$name <- nodePlot$name
      V(network)$color <- ifelse(V(network)$type == "metabolite", "#C6E3A1", "#FED7D9")
      
      pdf(paste("./Output/M4_MPNet/network_",job_id,"_time_",timen,".pdf",sep = ""), width = plot_size, height = plot_size)
      set.seed(1)
      plot.igraph(network, 
                  mark.col = rgb(0.8,0.8,0.8,0.25),
                  mark.border = "grey85",
                  vertex.lable =V(network)$name,
                  vertex.color = V(network)$color,
                  vertex.size = 30,
                  vertex.label.color = "black",
                  vertex.label.family = "sans",  
                  edge.color = "grey45", 
                  edge.weidth = linkPlot$edge_weight*1.5,
                  vertex.frame.color ="white",
                  margin = rep(0, times = 4) )
      dev.off()
    }else{
      # control the size of node
      daforsize <- data.frame("degree"=unique(nodePlot$degree),"size"=NA)
      if(max(daforsize$degree) <= 30 & min(daforsize$degree) >= 5){
        V(network)$size <- nodePlot$degree
      }else{
        step = (30-5)/max(c( nrow(daforsize)-1 , 1)) 
        daforsize <- daforsize[order(-daforsize$degree),]
        daforsize$size <- seq(30,5,-step)[1:nrow(daforsize)]
        nodePlot <- dplyr::left_join(nodePlot, daforsize,by="degree")
        V(network)$size <- nodePlot$size
        nodePlot <- nodePlot[,which(colnames(nodePlot) != "size")]
      }
      
      V(network)$name <- nodePlot$name
      V(network)$color <- ifelse(V(network)$type == "metabolite", "#C6E3A1", "#FED7D9")
      
      com <- list()
      for(i in unique(nodePlot$community)){
        com[[paste("com",i,sep = "")]] <- nodePlot[which(nodePlot$community == i),]$name
      }
      
      e <- igraph::as_edgelist(network,names=FALSE)
      l <- qgraph::qgraph.layout.fruchtermanreingold(e, vcount = vcount(network),
                                                     area = 6*(vcount(network)^2), repulse.rad = vcount(network)^3)
      
      
      pdf(paste("./Output/M4_MPNet/network_",job_id,"_time_",timen,".pdf",sep = ""), width = plot_size, height = plot_size)
      set.seed(1)
      plot.igraph(network, 
                  layout = l,
                  mark.groups = com,
                  mark.col = rgb(0.8,0.8,0.8,0.25),
                  mark.border = "grey85",
                  vertex.lable =V(network)$name,
                  vertex.color = V(network)$color,
                  vertex.size = V(network)$size,
                  vertex.label.color = "black",
                  vertex.label.family = "sans",  
                  edge.color = "grey45", 
                  edge.weidth = linkPlot$edge_weight*1.5,
                  vertex.frame.color ="white",
                  margin = rep(0, times = 4) )
      dev.off()
    }
    
  } 
   
}

# ------------------Test ----------------------------------------

 
set.seed(42)
proteins_test <- proID_ref_annotation[sample(1:5000,100 ,replace = F),2]
metabolites_test <- rio::import("./testda/testData2.txt")$x


res <- ProMetNetwork(proteins = proteins_test, metabolites = metabolites_test, 
                   protein_type = "Uniprot",job_id = "test42", timen = 0)


