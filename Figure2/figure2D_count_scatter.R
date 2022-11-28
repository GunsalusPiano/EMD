

########################
# R script for figure 2 d: Cell counts scatterplot
########################

# clear variables
rm(list=ls())

# load libraries
library(rmarkdown)
library(easypackages)
library(resample) 
libraries("assertthat","lattice","gplots","ggplot2","fields","readr","data.table","grid","MASS")  
libraries("plotly","gridExtra","Rmisc","Hmisc","doBy","caret","corrplot","GGally", "easypackages")  
library(reshape) 
library(plyr)
library(ggpubr)
library(tictoc)

   
  ########################
  # import the annotated well median data
  ########################

    # setwd("~/Directory.../.../")
    well_med_dat <- read_csv("raw_medians_fig2C.csv")
    well_med_dat = well_med_dat[,-1]

  #########################
  # graphing cell counts
  #########################

  names(well_med_dat)[12] = "cell_count"
  
  #########################
  # adding indices
  #########################
  
  well_med_dat$x = 1:1848
  
  #########################
  # control wells
  #########################
  
  ind_control = grep('control', well_med_dat$Treatment)
  annot.count.c = well_med_dat[ind_control,]
  
  #########################
  # summary of control well counts
  #########################
  
  control_range = range(annot.count.c$cell_count)
  control_range
  control_med = median(annot.count.c$cell_count)
  control_med
  count.r = range(well_med_dat$cell_count)

  #########################
  # check number of plates
  #########################
  
  names(well_med_dat)
  unique(well_med_dat$PlateNumber)

  #########################
  # remove empty treatments
  #########################
  
  unique(well_med_dat$Treatment)
  annot.count = well_med_dat[-grep("empty",well_med_dat$Treatment),]
  
  # adjusting the metadata - adding column for annotated dilution factor
  annot.count$Concentration = annot.count$Metadata_dil_factor
  unique(annot.count$Metadata_dil_factor)
  
  # "NULL"     "1"        "0.5"      "0.25"     "0.125"    "0.0625"   "0.03125"  "0.015625"
  unique(annot.count$Concentration)
  annot.count.p1 = annot.count
  annot.count.p1$Concentration[which(annot.count.p1$Concentration == "NULL")]="DMSO"
  annot.count.p1$Concentration[which(annot.count.p1$Concentration == "1")]="High"
  annot.count.p1$Concentration[which(annot.count.p1$Concentration == "0.5")]="High"
  annot.count.p1$Concentration[which(annot.count.p1$Concentration == "0.25")]="Med"
  annot.count.p1$Concentration[which(annot.count.p1$Concentration == "0.125")]="Med"
  annot.count.p1$Concentration[which(annot.count.p1$Concentration == "0.0625")]="Low"
  annot.count.p1$Concentration[which(annot.count.p1$Concentration == "0.03125")]="Low"
  annot.count.p1$Concentration[which(annot.count.p1$Concentration == "0.015625")]="Low"
  
    #########################
    # check the metadata
    #########################
  
    unique(annot.count.p1$Concentration)
    
    #########################
    # adding metadata
    #########################
  
    annot.count.p1$Treatment2=annot.count.p1$Treatment
    annot.count.p1$Treatment2[which(annot.count.p1$Treatment2 == "treated")]="Treatment wells"
    annot.count.p1$Treatment2[which(annot.count.p1$Treatment2 == "control")]="Control wells"
    annot.count.p1$Treatment1 = factor(annot.count.p1$Treatment2, levels = c("Treatment wells", "Control wells"))
  
  
    #########################
    # Counts plotting paramters
    #########################
    
    Y.text <- element_text( color = "black", size = 14, angle = 0)
    X.text <- element_text( color = "black", size = 16, angle = 0)
    a.text = element_text(size = 16, color = "black")
    t.text = element_text(size = 14, color = "black")
    strip.text.a = element_text(size = 20, color = "black")
    
    
    #########################
    # Cell counts plot showing both treatment and control
    #########################
    
      ggplot(annot.count.p1, aes(x, cell_count, colour = Concentration, shape = Treatment)) +
      geom_point(shape = 20, size = 2, show.legend = FALSE, alpha = .6) + 
      scale_color_manual(values=c("indianred4","steelblue4","lightblue","steelblue2")) + 
      scale_shape_manual(values=c(20, 20)) +
      ###
      geom_hline(yintercept= control_med, linetype ="dashed", color = "darkgrey", size = 1) +
      geom_hline(yintercept = control_range[1], linetype="dashed", color = "darkgrey") +
      geom_hline(yintercept = control_range[2], linetype="dashed", color = "darkgrey") +
      ###
      ylim(-10, 1200) +
      theme_classic() +
      theme(strip.text = strip.text.a,axis.title.x=element_blank(), axis.text.x = X.text, axis.text.y = Y.text, axis.title = a.text, title = t.text) + 
      labs(y = NULL, title = NULL) +
      scale_x_continuous(breaks = c(1, 309, 617, 925, 1233, 1541), 
                         labels=c("p1.r1", "p1.r2", "p1.r3", "p2.r1", "p2.r2", "p2.r3")) + 
      facet_wrap(~Treatment2)
    
    
  