

########################
# R script for figure  7 A: UMAP trajectories
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
library("scatterplot3d") 
library(rgl)
library(rafalib)

library(RColorBrewer)



########################
# Function for connecting treatments (compounds) by concentrations
########################


# set working directory

# setwd("~/Directory.../.../")
source("line_segments.R")


########################
# import umap dimensions (see fig 7A code for generating the umap dimensions)
########################

# setwd("~/Directory.../.../")
df <- read_csv("UMAP_dimensions_69feats_785treatments_new.csv")


#########################
# Preparing UMAP trajectory plots
#########################


########################
# Control samples 
########################

df.umap = df 

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
# A groupd of diverse compounds
########################

comp.set = c("digitoxin", "tolfenamic-acid", "bleomycin", "methotrexate", 
        "FCCP", "irinotecan", "nocodazole", "tanespimycin")


########################
mypalette1 <- brewer.pal(7, "Purples")
#image(1:7,1,as.matrix(1:7),col=mypalette1,xlab = "Purples (sequential)", ylab="", xaxt="n", yaxt="n", bty="n")
col.use.p = c(mypalette1[5], mypalette1[6]) 
########################

########################
mypalette2 <- brewer.pal(7, "Reds")
#image(1:7,1,as.matrix(1:7),col=mypalette2, xlab = "Reds (sequential)", ylab="", xaxt="n", yaxt="n", bty="n")
col.use.r = c(mypalette2[5], mypalette2[6]) 
########################

#######################################
mypalette3 <- brewer.pal(7, "Blues")
#image(1:7,1,as.matrix(1:7), col = mypalette3, xlab = "Blues (sequential)", ylab = "", xaxt = "n", yaxt = "n", bty="n")
col.use.b = c(mypalette3[5], mypalette3[6]) 
#######################################

########################
mypalette4 <- brewer.pal(7, "Greens")
#image(1:7,1, as.matrix(1:7), col = mypalette4, xlab = "Greens (sequential)", ylab = "", xaxt = "n", yaxt = "n", bty="n")
col.use.g = c(mypalette4[5], mypalette4[6]) 
########################


yes.trt = comp.set 


########################
# Compound name is yes.trt
Features_umap_df.drug <- Features_umap_df[which(Features_umap_df$Drug_name %in% yes.trt),]

########################
# umap.lines.select <- make_line_segments(Features_umap_df.drug)
########################
# dim1 and dim2 -- connecting lines
umap.lines.sub <- make_line_segments(Features_umap_df.drug)
# dim2 and dim3 -- connecting lines
umap.lines.sub1 <- make_line_segments1(Features_umap_df.drug)
########################




########################
# RED, GREEN, BLUE, PURPLE, RED2, GREEN2, BLUE2, PURPLE2
########################

 mycolors_highlight = c(col.use.p[1], "#117733", "#88CCEE", 
                        "#6699CC", col.use.p[2], "#882255", "#661100", "#999933")

# mycolors_highlight = c(col.use.p[1], col.use.g[1], col.use.b[1], 
#                        col.use.b[2], col.use.p[2], col.use.r[1], col.use.r[2], col.use.g[2])



########################
# UMAP Dim1 vs. Dim2
########################

ggplot() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  # theme(legend.position = "top") +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(size = 15)) +
  theme(axis.text.y = element_text(size = 15)) +
  theme(axis.title.x = element_text(size = 17)) +
  theme(axis.title.y = element_text(size = 17)) +
  labs(x = "\nUMAP D1", y = "UMAP D2") +
  # ylim(-5,3) + xlim(-5,3) +
  geom_vline( xintercept = 0, color = "grey") +
  geom_hline( yintercept = 0, color = "grey") +
  #
  geom_point(data = Features_umap_df.controls, mapping=aes(x=V1, y=V2), size= 2, shape = 17, colour = "gray") +
  geom_point(data = Features_umap_df, alpha = 0.8, stroke = 1.1, shape = 21, fill = "white", colour = "gray", mapping=aes(x=V1, y=V2,  size=(0.40*Metadata_dil_factor + 0.02))) +
  geom_curve(data = umap.lines, alpha = 0.2, colour = "gray", mapping = aes(x=x1, y=y1, xend=x2, yend=y2,  size=(0.65*Metadata_dil_factor)) ,curvature = -0.2 ) +
  # extra layer for colored paths for compoounds that effect cell cycle 
  geom_curve(data = umap.lines.sub, mapping=aes(x=x1, y=y1, xend=x2, yend=y2, colour = factor(Drug_name), size=(0.85*Metadata_dil_factor)), curvature = -0.2) +
  geom_point(data = Features_umap_df.drug, stroke = 1.15, shape = 21, fill = "white", mapping=aes(x=V1, y=V2, colour = factor(Drug_name), size=(0.40*Metadata_dil_factor + 0.02))) + 
  scale_color_manual(values = mycolors_highlight)




########################
# Dim2 vs. Dim3
########################

ggplot() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  # theme(legend.position = "top") +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(size = 15)) +
  theme(axis.text.y = element_text(size = 15)) +
  theme(axis.title.x = element_text(size = 17)) +
  theme(axis.title.y = element_text(size = 17)) +
  labs(x = "\nUMAP D2", y = "UMAP D3") +
  # ylim(-5,3) + xlim(-5,3) +
  geom_vline( xintercept = 0, color = "grey") +
  geom_hline( yintercept = 0, color = "grey") +
  
  # conntrols
  geom_point(data = Features_umap_df.controls, mapping=aes(x=V2, y=V3), size = 2, shape = 17, colour = "gray") +
  # points
  geom_point(data = Features_umap_df, alpha = 0.8, stroke = 1.1, shape = 21, fill = "white", colour = "gray", mapping=aes(x=V2, y=V3,  size=(0.40*Metadata_dil_factor + 0.02))) +
  # grey curves
  geom_curve(data = umap.lines1, alpha = 0.2, colour = "gray", mapping = aes(x=x1, y=y1, xend=x2, yend=y2,  size=(0.45*Metadata_dil_factor)) ,curvature = -0.2 ) +
  
  # extra layer for colored paths for compoounds that effect cell cycle 
  # dotted, longdash, dashed
  # size=(0.45*Metadata_dil_factor)
  
  geom_curve(data = umap.lines.sub1, alpha = 0.9, mapping=aes(x=x1, y=y1, xend=x2, yend=y2, colour = factor(Drug_name),  size=(0.45*Metadata_dil_factor)), curvature = -0.2) +
  geom_point(data = Features_umap_df.drug, stroke = 1.15, shape = 21, alpha = 0.8, fill = "white", mapping=aes(x=V2, y=V3, colour = factor(Drug_name), size=(0.40*Metadata_dil_factor + 0.02))) + 
  scale_color_manual(values = mycolors_highlight)




########################
# Legend
########################

col1=c("#117733", "#999933", col.use.p[1], col.use.p[2], "#88CCEE", "#6699CC", "#882255", "#661100")

plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')

labels = c("digitoxin", "tolfenamic-acid", "bleomycin", "methotrexate", 
           "FCCP", "irinotecan", "nocodazole", "tanespimycin")

legend("center", 
       legend = labels, 
       col = col1,
       pch = c(15, 15, 15, 15, 15, 15, 15, 15 ), 
       bty = "n", 
       pt.cex = 2, 
       cex = 1.2, 
       text.col = "black", 
       horiz = F , 
       inset = c(0.1, 0.1))

