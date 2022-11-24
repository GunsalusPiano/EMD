########################
# R script for figure  5 A: Merging replicates
########################

# clear variables
rm(list=ls())

########################
# Load libraries
########################
library(easypackages);library(resample);
libraries("assertthat","lattice","gplots","ggplot2","fields","readr","data.table","grid","MASS")  
libraries("plotly","gridExtra","Rmisc","Hmisc","doBy","caret","corrplot","GGally", "easypackages")  
library(reshape);
library(plyr);library(ggpubr);
library("scatterplot3d") ;library(rgl);library(rafalib);library(ggridges)


########################
# Import output single cell data for controls and vincristine treated cells  
########################



########################
# set working directory
# setwd("~/Documents/...")

dat.use = read_csv("cell_data_figure5a_vincristine_control.csv")


###
# grab the controls
###

ind.c = which(dat.use$Drug_name == "DMSO")

dat.c = dat.use[ind.c,]

###
# grab vincristine @ 20 um concentration
### 

ind.v = which(dat.use$drug_conc == "vincristine_20_1")

dat.vin = dat.use[ind.v,]



Y.text <- element_text( color = "black", size = 20, angle = 0)
X.text <- element_text( color = "black",  size = 20, angle = 0) 
a.text = element_text(size= 20, color = "black")
t.text = element_text(size= 10, color = "black")


#########
# Plot 1
#########

# The replicate treatments
a = ggplot(dat.vin, aes(x = ObjectArea_NUC_A, color = plate_well)) +
  
  geom_density(size = 1.2, alpha = .4) + 
  
  theme(axis.text.y = Y.text, axis.text.x = X.text, axis.title = a.text, title = t.text, plot.title = element_text(size=22) ) +
  
  theme_light() +
  # theme_classic() +
  
  labs(x = NULL, y = "Density", title = NULL) +
  
  theme(legend.position = 'none') + 
  
  scale_color_manual(values=c("darkviolet", "deeppink","darkturquoise"))

# Add on the control
a1 = a +  
  
  geom_density(data = dat.c, aes(x=ObjectArea_NUC_A , color = Treatment), size = 1.5, col = "black") 

a1


#########
# Plot 2
#########

# The replicate treatments
b = ggplot(dat.vin, aes(x = ObjectArea_NUC_A, color = Treatment)) +
  
  geom_density(size = 1.2, alpha = .4) + 
  
  theme(axis.text.y = Y.text, axis.text.x = X.text, axis.title = a.text, title = t.text, plot.title = element_text(size=22) ) +
  
  theme_light() +
  # theme_classic() +
  
  labs(x = NULL, y = "Density", title = NULL) +
  
  theme(legend.position = 'none') +
  scale_color_manual(values=c('treated' = "#E69F00",'control' = "black")) 



# Add on the control
b1 = b +  
  
  geom_density(data = dat.c, aes(x=ObjectArea_NUC_A , color = Treatment), size = 1.5, col = "black") 

b1




dat.all = rbind(dat.vin,dat.c)


#########
# Plot 3
#########

# CDF of merged treatments
cdf = ggplot(dat.all, aes(x = ObjectArea_NUC_A, color = Treatment)) + 
  stat_ecdf(size = 1.2) +
  
  theme_light() +
  theme(legend.position = 'none') +
  theme(axis.text.y = Y.text, axis.text.x = X.text, axis.title = a.text, 
        title = t.text, plot.title = element_text(size=22) ) +
  
  labs(x = NULL, y = "ECDF", title = NULL) +
  scale_color_manual(values=c('treated' = "#E69F00",'control' = "black")) 



grid.arrange(a1,b1,cdf, nrow=1)
