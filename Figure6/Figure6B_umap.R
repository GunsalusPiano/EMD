

########################
# R script for figure  6 B: UMAP 
########################


########################
# Clear variables
########################

rm(list=ls())


########################
# Load libraries
########################

library(easypackages)
library(resample) 
library(caroline)
libraries("assertthat","lattice","gplots","ggplot2","fields","readr","data.table","grid","MASS")  
libraries("plotly","data.table","gridExtra","Rmisc","Hmisc","doBy","caret","corrplot","GGally", "easypackages")  
library(reshape) 
library(plyr)
library(ggpubr)
library(Rtsne)
library(tsne)
library("scatterplot3d") 
library(rgl)
library(rafalib)
library(tictoc)
library(ggridges)
library(plyr)
library(RColorBrewer)
library(umap)


########################
# Function for connecting treatments (compounds) by concentrations
########################

# setwd("~/Directory.../.../")
source("line_segments.R")


########################
# Import EMD profile
########################

# setwd("~/Directory.../.../")
df.emd.raw  = read_csv("fullfeature_EMDprofile_785treatments.csv")


########################
# Reduce profile to optimal features
########################

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


dfraw = df.emd.raw 
df69 = dfraw[,which(colnames(dfraw) %in% f.keep)]
df.a = dfraw[, 2:17] 
df69.all = cbind(df.a, df69)

names(df69.all)


#########################
# Find first feature
#########################

df.full = df69.all 
f_start_col = which("ObjectShapeP2A_NUC_A"   ==  names(df.full))
f_start_col

df.full.num = df.full[,f_start_col:ncol(df.full)]
df.num = df.full.num

#########################
# Log transform the feature data 
#########################

# natural log
dlog = log(df.num) 

# min-max scale the features to [0,1]
max.log <- max(dlog, na.rm = TRUE)
max.log

min.log <- min(dlog, na.rm = TRUE)
min.log

dlog.sc = (dlog - min.log)/(max.log - min.log)

rownames(dlog.sc) = df.full$drug_conc

#########################
# Full feature set minmax_scaled(Log_scaled)
d.sc = cbind(df.full[,1:16], dlog.sc)
d.sc.temp = d.sc
#########################

#########################
# create data frame d.new
# convert profile to matrix
d.new = as.matrix(dlog.sc)
rownames(d.new) = d.sc$drug_conc


#########################
# run the umap
#########################

custom.config = umap.defaults
custom.config$n_components <- 3

Features_umap <- umap(d.new,config = custom.config)
dim(Features_umap$layout)

Features_umap_df <- as.data.frame(Features_umap$layout)

# adding annotations 
Features_umap_df1 <- cbind(Features_umap_df,d.sc[,c("Drug_name", "Metadata_dil_factor", "Metadata_concentration_uM", "cell_count_POC_A", "drug_conc")])
                                                     
#########################
# Graph using Drug name
plot_ly(Features_umap_df1,x = ~V1, y= ~V2, z = ~V3, type = "scatter3d", color = ~Drug_name, mode = "markers")


#########################
# UMAP trajectory plots prep
#########################


########################
# Control samples 
########################

df.umap = Features_umap_df1

Features_umap_df.controls = df.umap[grep("DMSO", df.umap$Drug_name, invert = FALSE),]

########################
# Treatment samples
########################

Features_umap_df.t <- df.umap[grep("DMSO", df.umap$Drug_name, invert = TRUE),]

########################
# dim1 and dim2 -- connecting lines
########################

umap.lines <- make_line_segments(Features_umap_df.t)

########################
# dim2 and dim3 -- connecting lines
########################

umap.lines1 <- make_line_segments1(Features_umap_df.t)

########################
# treatments umap dimensions
########################

Features_umap_df = Features_umap_df.t

########################
# adding dilution factor meta as numeric data to the line segment data
# dim1 dim2
########################

str(umap.lines)
umap.lines$Metadata_dil_factor = as.numeric(as.character(umap.lines$Metadata_dil_factor))
str(umap.lines)
unique(umap.lines$Metadata_dil_factor)

########################
# adding dilution factor meta as numeric data to the line segment data
# dim2 dim3
########################

str(umap.lines1)
umap.lines1$Metadata_dil_factor = as.numeric(as.character(umap.lines1$Metadata_dil_factor))
str(umap.lines1)
unique(umap.lines1$Metadata_dil_factor)

str(Features_umap_df)
Features_umap_df$Metadata_dil_factor = as.numeric(as.character(Features_umap_df$Metadata_dil_factor))
str(Features_umap_df)
unique(umap.lines$Metadata_dil_factor)


########################
# UMAP plot
# Dim1 vs. Dim2
########################

ggplot() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(size = 15)) +
  theme(axis.text.y = element_text(size = 15)) +
  theme(axis.title.x = element_text(size = 17)) +
  theme(axis.title.y = element_text(size = 17)) +
  labs(x = "\nUMAP D1", y = "UMAP D2") +
  
  geom_vline( xintercept = 0, color = "grey") +
  geom_hline( yintercept = 0, color = "grey") +
  geom_point(data = Features_umap_df.controls, mapping = aes(x = V1, y = V2, colour = cell_count_POC_A), size = 2, shape = 17) +
  geom_point(data = Features_umap_df, alpha = 0.8, stroke = 1.1, shape = 21, fill = "white", mapping = aes(x = V1, y = V2, colour = cell_count_POC_A, size=(0.40*Metadata_dil_factor + 0.02))) +
  scale_color_gradient(low = "red", high = "blue") +
  geom_curve(data = umap.lines, alpha = 0.2, colour = "gray", mapping = aes(x = x1, y = y1, xend = x2, yend = y2,  size = (0.65*Metadata_dil_factor)) ,curvature = -0.2 )





########################
# UMAP plot
# Dim2 vs. Dim3
########################

ggplot() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(size = 15)) +
  theme(axis.text.y = element_text(size = 15)) +
  theme(axis.title.x = element_text(size = 17)) +
  theme(axis.title.y = element_text(size = 17)) +
  labs(x = "\nUMAP D2", y = "UMAP D3") +
  
  geom_vline( xintercept = 0, color = "grey") +
  geom_hline( yintercept = 0, color = "grey") +
  geom_point(data = Features_umap_df.controls, mapping=aes(x=V2, y=V3, colour = cell_count_POC_A), size= 2, shape = 17) +
  geom_point(data = Features_umap_df, alpha = 0.8, stroke = 1.1, shape = 21, fill = "white", mapping=aes(x=V2, y=V3, colour = cell_count_POC_A, size=(0.40*Metadata_dil_factor + 0.02))) +
  scale_color_gradient(low ="red", high ="blue") +
  geom_curve(data = umap.lines1, alpha = 0.2, colour = "gray", mapping = aes(x=x1, y=y1, xend=x2, yend=y2,  size=(0.65*Metadata_dil_factor)) ,curvature = -0.2 )










