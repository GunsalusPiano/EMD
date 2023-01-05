########################
# R script for figure 3D: Pre and post corrected data
########################


  rm(list=ls())

########################
# Load libraries
########################

  library(easypackages)
  library(resample) 
  library(caroline)
  libraries("assertthat","lattice","gplots","ggplot2","fields","readr","data.table","grid","MASS")  
  libraries("plotly","gridExtra","Rmisc","Hmisc","doBy","caret","corrplot","GGally", "easypackages")  
  library(reshape) 
  library(plyr)
  library(ggpubr)
  library(RColorBrewer)


    setwd("~/Documents/NYU_work/Sept_17_2018/R_code_b/Manuscript_Comm_Bio_scripts/Figure3")

    # setwd("~/Directory.../.../")
    df.new <- read_csv("cell_data_fig3D.csv")

    df = df.new
    
    ########################
    # ggplot parameters 
    ########################
    Y.text <- element_text(color = "black", size = 15, angle = 0)
    X.text <- element_text(color = "black", size = 15, angle = 0)
    a.text = element_text(size = 15, color = "black")
    t.text.o = element_text(size = 15, color = "red")
    t.text.n = element_text(size = 16 ,color = "green")
    strip.text.a = element_text(size = 14, color = "purple")

    
    
    ########################
    # raw cell feature distributions
    ########################

    a2 = 
      ggplot(df, aes(x=FeatureX_raw_cells, color=PlateNumber)) +
     geom_density(size = 1.1) +
     scale_color_brewer(palette="Dark2") + 
     theme_classic() + 
     theme(legend.position="none") +
     labs(x = "Total Nucleus Intensity \n (raw data)",y = "Per plate density") +
     theme(axis.text.x = X.text, 
        axis.text.y = Y.text, axis.title = a.text, title = t.text.o)


    ########################
    # adjusted cell feature distributions
    ########################
    

    b2 = 
      ggplot(df, aes(x=FeatureX_adj_cells, color=PlateNumber)) +
    geom_density(size = 1.1) +
    scale_color_brewer(palette="Dark2") + 
    theme_classic() + 
    theme(legend.position="none") +
    labs(x = "Total Nucleus Intensity \n (adjusted for positional effects)", y = NULL) +
    theme(axis.text.x = X.text, axis.text.y = Y.text, axis.title = a.text, title = t.text.o)
        

    ########################
    # standardized cell feature distributions
    ########################
    

    c2 = 
      ggplot(df, aes(x=FeatureX_norm_cells, color=PlateNumber)) +
     geom_density(size = 1.1) +
     scale_color_brewer(palette="Dark2") + 
     theme_classic() +  
     theme(legend.position = "none") +
     theme(legend.title = element_text(colour="black", size=14)) +
     theme(legend.text = element_text(colour="black", size=14)) +
     labs(x = "Total Nucleus Intensity \n (adjusted for plate effects)", y = NULL) +
     theme(axis.text.x = X.text, axis.text.y = Y.text, axis.title = a.text, title = t.text.o)


    ########################
    # Combined figures for Figure 3D
    ########################
    
  
    grid.arrange(a2,b2,c2,nrow=1)


