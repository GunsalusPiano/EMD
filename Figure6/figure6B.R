
########################
# R script for figure  6 B: Hierarchical clustering and dimension reduction.
########################

########################
# Clear variables
# Set working directory
########################

rm(list=ls())

# setwd("~/Directory.../.../")


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
library("scatterplot3d") 
library(rgl)
library(rafalib)
library(tictoc)




########################
# Import annotated UMAP dimensions
########################

setwd("/Users/pearsy02/Documents/NYU_work/Sept_17_2018/R_code_b/Manuscript_Comm_Bio_scripts")

########################
# Function for connecting treatments (compounds) by concentrations
########################

source("line_segments.R")

########################
# Import annotated UMAP dimensions
########################


# df.umap = read_csv("~/Directory/.../umap_dims.csv")


########################
# Control samples 
########################
Features_umap_df.controls = df.umap[grep("DMSO", df.umap$Drug_name, invert = FALSE),]

########################
# Treatment samples
########################
temp.t <- df.umap[grep("DMSO", df.umap$Drug_name, invert = TRUE),]
Features_umap_df.t <- temp.t

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
  theme(legend.position = "top") +
  theme(axis.text.x = element_text(size = 15)) +
  theme(axis.text.y = element_text(size = 15)) +
  theme(axis.title.x = element_text(size = 17)) +
  theme(axis.title.y = element_text(size = 17)) +
  labs(x = "\nUMAP D1", y = "UMAP D2") +
  
  geom_vline( xintercept = 0, color = "grey") +
  geom_hline( yintercept = 0, color = "grey") +
  geom_point(data = Features_umap_df.controls, mapping = aes(x = V1, y = V2, colour = Cell_count_poc), size = 2, shape = 17) +
  geom_point(data = Features_umap_df, alpha = 0.8, stroke = 1.1, shape = 21, fill = "white", mapping = aes(x = V1, y = V2, colour = Cell_count_poc, size=(0.40*Metadata_dil_factor + 0.02))) +
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
  geom_point(data = Features_umap_df.controls, mapping=aes(x=V2, y=V3, colour = Cell_count_poc), size= 2, shape = 17) +
  geom_point(data = Features_umap_df, alpha = 0.8, stroke = 1.1, shape = 21, fill = "white", mapping=aes(x=V2, y=V3, colour = Cell_count_poc, size=(0.40*Metadata_dil_factor + 0.02))) +
  scale_color_gradient(low ="red", high ="blue") +
  geom_curve(data = umap.lines1, alpha = 0.2, colour = "gray", mapping = aes(x=x1, y=y1, xend=x2, yend=y2,  size=(0.65*Metadata_dil_factor)) ,curvature = -0.2 )
