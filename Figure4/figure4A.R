


########################
# R script for figure  4 A: Replicate distributions - Control and Treated cells
########################




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
library("FactoClass")
library(rafalib)
library(tictoc)
library(ggridges)
library(plyr)



  ########################
  # Import normalized single cell data  
  ########################


  
setwd("/Users/pearsy02/Documents/NYU_work/Sept_17_2018/R_code_b/Manuscript_Comm_Bio_scripts/Data_files/Fig4")

########################
# set working directory
# setwd("~/Documents/...")
  

  ########################
  # import cell data 
  dat.c <- read_csv("cell_data_figure4A_control.csv")
  dat.t1 <- read_csv("cell_data_figure4A_treatment1.csv")
  dat.t2 <- read_csv("cell_data_figure4A_treatment2.csv")
  

########################  
# check metadata
  
unique(dat.c$Plate) # "Plate1_drugs" "Plate2_drugs"
unique(dat.c$Plate_number)
unique(dat.c$PlateNumber) # "1A1" "1A3" "1A2" "2A1" "2A2" "2A3"

########################
# renaming the dataframe
new_ctrl = dat.c

########################
# add meta data for the DMSO cells
new_ctrl$PlateRep = "Replicate1"
ind1 = which(new_ctrl$PlateNumber == "1A1")
new_ctrl$PlateRep[ind1] = "Replicate1"
ind2 = which(new_ctrl$PlateNumber == "1A2")
new_ctrl$PlateRep[ind2] = "Replicate2"
ind3 = which(new_ctrl$PlateNumber == "1A3")
new_ctrl$PlateRep[ind3] = "Replicate3"

ind4 = which(new_ctrl$PlateNumber == "2A1")
new_ctrl$PlateRep[ind4] = "Replicate1"
ind5 = which(new_ctrl$PlateNumber == "2A2")
new_ctrl$PlateRep[ind5] = "Replicate2"
ind6 = which(new_ctrl$PlateNumber == "2A3")
new_ctrl$PlateRep[ind6] = "Replicate3"

unique(new_ctrl$PlateRep)
unique(new_ctrl$plate_well)

new_ctrl$plate_replicate = paste(new_ctrl$Plate, new_ctrl$PlateRep, sep = "_")
unique(new_ctrl$plate_replicate)


########################
# add meta data for the chemically treated cells
########################

# renaming the dataframe
new_trt = rbind(dat.t1, dat.t2)

new_trt$PlateRep = "Replicate1"

ind1 = which(new_trt$PlateNumber == "1A1")
new_trt$PlateRep[ind1] = "Replicate1"
ind2 = which(new_trt$PlateNumber == "1A2")
new_trt$PlateRep[ind2] = "Replicate2"
ind3 = which(new_trt$PlateNumber == "1A3")
new_trt$PlateRep[ind3] = "Replicate3"

ind4 = which(new_trt$PlateNumber == "2A1")
new_trt$PlateRep[ind4] = "Replicate1"
ind5 = which(new_trt$PlateNumber == "2A2")
new_trt$PlateRep[ind5] = "Replicate2"
ind6 = which(new_trt$PlateNumber == "2A3")
new_trt$PlateRep[ind6] = "Replicate3"

unique(new_trt$PlateRep)
unique(new_trt$plate_well)

########################
# add meta data to differentiate cells from replicate plates
########################

new_trt$plate_replicate = paste(new_trt$Plate, new_trt$PlateRep, sep = "_")
unique(new_trt$plate_replicate)

########################
# ggplot parameters
########################

Y.text <- element_text( color = "black", size =10, angle = 0)
X.text <- element_text( color = "black", size = 10, angle = 0)
a.text = element_text(size=15, color = "indianred4")
t.text = element_text(size=13 ,color = "indianred4")
strip.text.a = element_text(size = 17, color = "indianred4")

########################
# Plot the DMSO distributions
########################

  ggplot(new_ctrl, aes(x = ObjectTotalInten_NUC_A, color = WellId, group = WellId)) +
  geom_density(size = .2, alpha = .4) +
  scale_colour_grey() +
  theme_bw()+
  scale_y_continuous(limits=c(-0.05,1.5), breaks=c(0, 0.25, 0.5,
                                                   0.75, 1, 1.25, 1.5), expand = c(0,0)) +
  scale_x_continuous(limits=c(-5,12), breaks=c(-4,-2,0, 2, 4,
                                               6, 8, 10, 12), expand = c(0,0)) +
  theme(strip.background = element_rect(fill="white"), strip.text = strip.text.a, 
        axis.text.x = X.text, axis.text.y = Y.text, axis.title = a.text, title = t.text) +
  labs(x = "Total Nucleus Intensity", y = NULL, title = NULL) +
  theme(legend.position = 'none') +  
  facet_wrap(~PlateRep, nrow = 3)


########################
# ggplot parameters for the treated cells
########################

strip.text.a = element_text(size = 17, color = "steelblue")
a.text = element_text(size=15, color = "steelblue")

########################
# Plot the Treatment distributions
########################

  ggplot(new_trt, aes(x = ObjectTotalInten_NUC_A, color = plate_well, group = plate_well)) +
  geom_density(size = .2, alpha = .4)+
  scale_colour_grey() +
  theme_bw()+
  scale_y_continuous(limits=c(-0.05,1.5), breaks=c(0, 0.25, 0.5,
                                                   0.75, 1, 1.25, 1.5), expand = c(0,0)) +
  scale_x_continuous(limits=c(-5,12), breaks=c(-4,-2,0, 2, 4,
                                               6, 8, 10, 12), expand = c(0,0)) +
  theme(strip.background = element_rect(fill="white"), strip.text = strip.text.a, 
        axis.text.x = X.text, axis.text.y = Y.text, axis.title = a.text, title = t.text) +
  labs(x = "Total Nucleus Intensity", y = NULL, title = NULL) +
  theme(legend.position = 'none') +  facet_wrap(~PlateRep, nrow=3)