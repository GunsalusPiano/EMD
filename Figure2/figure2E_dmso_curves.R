

########################
# R script for figure 2e: Per well density of DMSO treated cells
########################


# clear variables
rm(list=ls())


########################
# load libraries
library(easypackages)
library(ggridges)
libraries("rmarkdown","resample","assertthat","lattice","gplots","ggplot2",
          "fields","readr","data.table","grid","MASS", "reshape", "plyr", "ggpubr", "tictoc")  
libraries("plotly","gridExtra","Rmisc","Hmisc","doBy","caret","corrplot","GGally", "easypackages")  


########################
# import single cell data


# setwd("~/Directory.../.../")

cell.dat <- read_csv("cell_data_DMSO_figure2E.csv")



######################
# creating new annotation column called "plate_well"
cell.dat$plate_well = paste(cell.dat$PlateNumber, cell.dat$WellId, sep = "_")
unique(cell.dat$plate_well)
######################


  Y.text <- element_text( color = "black", size = 12, angle = 0)
  X.text <- element_text( color = "black", size = 12, angle = 0)
  a.text = element_text(size = 15, color = "black")
  t.text = element_text(size = 15, color = "black")



  p.c = 
  ggplot(cell.dat, aes(ObjectTotalIntenCh1)) +
  geom_density(color = "indianred3", size = 0.2, alpha = 0.3) + 
  geom_density(data = cell.dat, aes(ObjectTotalIntenCh1, group = plate_well), color = "indianred3", size = 0.2, alpha = 0.6) + 
  geom_density(data = cell.dat, aes(ObjectTotalIntenCh1), color = "black", size = 0.9, alpha = 1) + 
  
  theme(legend.position ='none') + 
  # scale_color_grey() +
  theme_classic() +
  theme(axis.text.y = Y.text, axis.text.x = X.text, axis.title = a.text,
        title = t.text, plot.title = element_text(size=14) ) +
  labs(x = "Total Nucleus Intensity", title = "Control samples", y = NULL)


  p.c +
  scale_y_continuous(limits=c(0, .00000020), breaks=c(0, .00000005,
                                                      .00000010, .00000015, .00000020), expand = c(0,0)) +
  scale_x_continuous(limits=c(0, 55297885), breaks=c(0,10000000,20000000, 30000000, 40000000, 50000000,60000000), expand = c(0,0)) 


