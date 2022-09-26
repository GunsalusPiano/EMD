

  ########################
  # R script for figure  6A: Hierarchical clustering and dimension reduction.
  ########################



  ########################
  # Clear variables
  # Set working directory
  ########################

  rm(list=ls())

  # setwd("~/Directory.../.../")


  ######################
  # load the libraries
  ######################
  
  library(easypackages)
  library(resample) 
  library(caroline)
  libraries("assertthat","lattice","gplots","ggplot2","fields","readr","data.table","grid","MASS")  
  libraries("plotly","data.table","gridExtra","Rmisc","Hmisc","doBy","caret","corrplot","GGally", "easypackages")  
  library(reshape) 
  library(plyr)
  library(ggpubr)
  library(Rtsne)
  library("scatterplot3d") 
  library(rafalib)
  library(rgl)
  library(rafalib)
  library(tictoc)
  library(RColorBrewer)

  



  ######################
  # Import the raw features - full set, unscaled
  ######################


  # df.emd = read_csv("~/Directory/.../EMD_rawprofile.csv")

  ######################
  # annotations
  names(df.emd[2:23])
  
  ######################
  # treatments
  unique(df.emd$drug_conc)

  ######################
  # number of features is 174 raw
  length(df.emd[24:197]) 
  
  ######################
  # number of control-DMSO rows
  length(which(df.emd$Metadata_type == "control"))
  
  ######################
  # number of treated rows
  length(which(df.emd$Metadata_type == "treated"))

 
  ######################
  # Identify most active features 
  ######################
  


f.keep = c("ObjectShapeP2A_NUC_A",       "ObjectShapeLWR_NUC_A",       "ObjectTotalInten_NUC_A",     "ObjectAvgInten_NUC_A",      
           "RingTotalInten_RNA",         "CircSpotTotalInten_RNA",     "CircSpotTotalArea_RNA",      "CircSpotAvgArea_RNA",       
           "CircSpotCount_RNA",          "TotalInten_RNA",             "RingTotalInten_PMG",         "RingAvgInten_PMG",          
           "RingSpotTotalInten_PMG",     "RingSpotAvgInten_PMG",       "RingSpotAvgArea_PMG",        "RingSpotCount_PMG",        
           "CircSpotTotalArea_PMG",      "CircSpotAvgArea_PMG",        "CircSpotCount_PMG",          "RingTotalInten_MITO",      
           "RingSpotTotalInten_MITO",    "RingSpotAvgArea_MITO",       "CircSpotTotalInten_NUC_Tex", "CircSpotTotalArea_NUC_Tex", 
           "CircSpotAvgArea_NUC_Tex",    "CircSpotCount_NUC_Tex",      "ObjectShapeP2A_NUC_B",       "ObjectShapeLWR_NUC_B",      
           "ObjectTotalInten_NUC_B",     "ObjectSize_NUC_B",           "RingTotalInten_Perox",       "RingAvgInten_Perox",        
           "RingSpotTotalInten_Perox",   "RingSpotTotalArea_Perox",    "RingSpotCount_Perox",        "TotalInten_Perox",          
           "AvgInten_Perox",             "RingTotalInten_Lipids",      "RingAvgInten_Lipids",        "RingSpotAvgInten_Lipids",   
           "RingSpotAvgArea_Lipids",     "RingSpotCount_Lipids",       "TotalInten_Lipids",          "AvgInten_Lipids",           
           "ObjectArea_NUC_C1",          "ObjectShapeLWR_NUC_C1",      "RingTotalInten_ER" ,         "RingSpotTotalInten_ER",     
           "RingSpotAvgInten_ER",        "RingSpotCount_ER",           "CircSpotAvgInten_ER",        "CircSpotTotalArea_ER",      
           "CircSpotAvgArea_ER",         "CircSpotCount_ER",           "TotalInten_ER",              "AvgInten_ER",               
           "ObjectShapeP2A_NUC_C2",      "ObjectShapeLWR_NUC_C2",      "ObjectSize_NUC_C2" ,         "RingTotalInten_TUB",        
           "RingSpotAvgInten_TUB",       "RingSpotTotalArea_TUB",      "RingSpotAvgArea_TUB",        "RingSpotCount_TUB",         
           "RingTotalInten_ACT",         "RingSpotTotalInten_ACT",     "RingSpotTotalArea_ACT",      "RingSpotCount_ACT",         
           "TotalInten_ACT" )


  dfraw = df.emd
  df69 = dfraw[,which(colnames(dfraw) %in% f.keep)]
  # annotations
  df.a = dfraw[,1:23] 
  # 69 features
  df69.all = cbind(df.a, df69)

  ######################
  # find first feature, column 50
  f_start_col = which("ObjectShapeP2A_NUC_A"   ==  names(df69.all))
  f_start_col

  ######################
  # numeric dataframe
  df69.num = df69.all[,f_start_col:ncol(df69.all)]

  ######################
  # Log transform the feature data again
  ######################
  
  dlog = log(df69.num)
  
  ######################
  # Scale the data to the range [0,1]
  ######################
  
  max.log <- max(dlog, na.rm=TRUE)
  max.log

  min.log <- min(dlog, na.rm=TRUE)
  min.log

  dlog.sc = (dlog-min.log)/(max.log - min.log)
  
  ######################
  # Put annotations back
  ######################
  
  rownames(dlog.sc) = df69.all$drug_conc
  d.sc = cbind(df69.all[,1:23], dlog.sc)
  d.sc.temp = d.sc


  ######################
  # create data frame d.new
  # convert profile to matrix
  
  d.new = as.matrix(dlog.sc)
  rownames(d.new) = df69.all$drug_conc

  

  ######################
  # Coloring of heatmap features (columns) and rows (treatments)
  ######################
  
  
  
  col.one = brewer.pal(11, "BrBG")

  channel_colors <- unlist(lapply(colnames(d.new),function(x){
  
  if(grepl("NUC_A",x))  col.one[1]        
  else if (grepl("RNA",x)) col.one[3]     
  else if (grepl("PMG",x))  col.one[4]    
  else if (grepl("MITO",x)) col.one[5]    
  else if (grepl("NUC_Tex",x)) col.one[1] 
  
  else if (grepl("NUC_B",x)) col.one[1]   
  else if (grepl("Lyso",x)) col.one[6]    
  else if (grepl("Perox",x)) col.one[7]   
  else if (grepl("Lipids",x)) col.one[8]  
  
  else if (grepl("NUC_C1",x)) col.one[2]  
  else if (grepl("ER",x)) col.one[9]     
  
  else if (grepl("NUC_C2",x)) col.one[1]  
  else if (grepl("TUB",x))    col.one[10] 
  else if (grepl("ACT",x))    col.one[11] 
  else 'lightgrey'
  
}))

  
  ######################
  # Coloring of heatmap rows (treatments)
  # Controls samples are yellow, and treatment samples blue.
  ######################
  
  treatment_colors <- unlist(lapply(df69.all$Treatment,function(x){
  
  if(grepl("control",x))   'yellow' 
  else if (grepl("treated",x))  'royalblue' 
  else 'lightgrey'
  
  }))



  ######################
  # Hierarchical clustering
  # Visualize with heatmap
  ######################
  
  
heatmap.2(d.new, Colv = T, Rowv = T, trace = "none", density ="none", col = redgreen, cexRow=.2, cexCol=0.3, 
          ColSideColors = channel_colors, RowSideColors = treatment_colors,
          hclust=function(x) hclust(x,method="average"))


