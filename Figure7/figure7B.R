
########################
# R script for figure  7 B: Phenotypic characterization
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
# Import cell counts data
########################

# df = read_csv("~/Directory/.../cell_counts.csv")

########################
# Control and treatment dataframes
########################

check1.c = df[grep("control", df$Treatment, invert = FALSE),]        
check1.t = df[grep("control", df$Treatment, invert = TRUE),]    


########################
# Summary of controls
########################

control.numeric = check1.c$cell_counts

control.mean = mean(control.numeric) 
control.mean

control.min = min(control.numeric) 
control.min

control.max = max(control.numeric)
control.max


########################
# Pick the subset of compounds
########################

sub = c("tolfenamic-acid", "methotrexate", "irinotecan", "tanespimycin")
drug.sub = which(check1.t$Drug_name %in% sub)

########################
# new data frame for subset of compounds
count.df = check1.t[drug.sub,]

########################
# update dataframe for ordering compounds in the plot
count.df$Drug_name_order <- factor(count.df$Drug_name, levels = sub)
count.df$label_num = as.numeric(count.df$label)


########################
mypalette1 <- brewer.pal(7, "Purples")
image(1:7,1,as.matrix(1:7),col=mypalette1,xlab = "Purples (sequential)", ylab="", xaxt="n", yaxt="n", bty="n")
col.use.p = c(mypalette1[5], mypalette1[6]) 
########################

########################
mypalette2 <- brewer.pal(7, "Reds")
image(1:7,1,as.matrix(1:7),col=mypalette2, xlab = "Reds (sequential)", ylab="", xaxt="n", yaxt="n", bty="n")
col.use.r = c(mypalette2[5], mypalette2[6]) 
########################

########################
mypalette3 <- brewer.pal(7, "Blues")
image(1:7,1,as.matrix(1:7), col = mypalette3, xlab = "Blues (sequential)", ylab = "", xaxt = "n", yaxt = "n", bty="n")
col.use.b = c(mypalette3[5], mypalette3[6]) 
########################

########################
mypalette4 <- brewer.pal(7, "Greens")
image(1:7,1, as.matrix(1:7), col = mypalette4, xlab = "Greens (sequential)", ylab = "", xaxt = "n", yaxt = "n", bty="n")
col.use.g = c(mypalette4[5], mypalette4[6]) 
########################


########################
# Colors for count dose response plot
########################

col.all = c(mypalette4[5], mypalette4[6], mypalette1[6], mypalette1[7], mypalette3[5], mypalette3[3],   mypalette2[5], mypalette2[6])


########################
# Plotting parameters
########################

Y.text <- element_text( color = "black", size =18, angle = 0)
X.text <- element_text( color = "black", size = 14, angle = 0)
a.text = element_text(size = 18, color = "black")
t.text = element_text(size = 18, color = "black")
strip.text.a = element_text(size = 24, color = "black")



p1 = 
  
  ggplot(count.df) + 
  aes(x = label_num, y = cell_counts, group = PlateNumber) +
  geom_line(linetype = 1, aes(colour = Drug_name), size = .8) +

  geom_point(aes(color= Drug_name), size=1.5, shape =  15) + ylim(0,1250) +
  scale_color_manual(values=c( "tolfenamic-acid" = col.all[2], "methotrexate" = col.all[4], 
                               "irinotecan" = col.all[5], "tanespimycin" = col.all[8])) +

  
  labs(y = "Raw Cell Count", x =  "Dilution Factor") +
   
  theme_light() +

  geom_hline(yintercept = control.min, linetype = "dashed", size = .5) +
  geom_hline(yintercept = control.mean, linetype = "dashed", size = 1.5) +
  geom_hline(yintercept = control.max, linetype = "dashed", size = .5) +
  
  theme(legend.position='none') +
  scale_x_continuous(breaks = c(1, 2, 3, 4, 5, 6, 7), 
                     labels=c("1/64", "1/32", "1/16", "1/8", "1/4", "1/2", "1")) +
  theme(strip.text = strip.text.a, axis.text.x = X.text, 
        axis.text.y = Y.text, axis.title = a.text, title = t.text) +
  theme(strip.background = element_rect(colour = "grey", fill = "white", size = 2, linetype = 1)) +
  facet_wrap(~ Drug_name_order, nrow = 1)




print(p1)




