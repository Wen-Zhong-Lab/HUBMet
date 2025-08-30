# example for visualization in MSEA

source("./Script/visualization.R") # functions for visualization

testda <- rio::import( "./testda/testda_hubmet_0709.txt" )
testda2 <- rio::import( "./testda/testda_public_0709.txt" )
 

hubres <- rio::import( "./Output/M2_MSEA/HUBMet_Disease_MSEA_forRvisu.txt" )
  
barplot_msea(hubres, "Disease")
  
msea_heatmap(hubres,hmdb_list = testda$metID,hmdb_list_value = testda$fc,database = "Disease")
  
msea_plot_ES_bar(database = "Disease",hubres$Term[1],hmdb_list = testda$metID,hmdb_list_value = testda$fc,res_msea = hubres)
 

pathview_kegg_smpdb_reactome("smpdb",pathway_ID = "SMP0000188", 
                             hmdb_list = testda2$metID,
                             hmdb_list_value = testda2$fc )

 