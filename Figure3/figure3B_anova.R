
########################
# R script for figure 3B: Two way anova analysis
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







# setwd("~/Directory.../.../")

dat <- read_csv("anova_output_fig3b.csv")



########################
# color code intensity and non-intensity features
########################

rr.rc = range(dat$P_c_log, dat$P_r_log)

low = rr.rc[1] - 3 
high = rr.rc[2] + 3


########################
# ggplot parameters
########################

Y.text <- element_text( color = "black", size = 12, angle = 0)
X.text <- element_text( color = "black", size = 12, angle = 0)
a.text = element_text(size = 18, color = "black")
t.text = element_text(size = 16 ,color = "black")
strip.text.a = element_text(size = 18, color = "black")


########################
# Figure 3B Summary of two-way ANOVA test for row effects
########################
  
  ggplot(dat) + 
  aes(x = plate_rep1, y = log_neg, group = feature, colour = feature_type) +
  geom_line(linetype = 1, size = .6) + 
  geom_hline(yintercept=-log(0.00001), linetype = "dashed", color = "black", size = 0.6) +
  geom_hline(yintercept=-log(0.001), linetype = "dotdash", color = "black", size = 0.6) +
  
  ylim(-high,-low) +
  theme_bw() +
  theme(legend.position = "none") +
  
  labs(x = NULL, y = "-log(p-value)") + 
  
  theme(strip.background = element_rect(fill="white"), strip.text = strip.text.a, axis.text.x = X.text, axis.text.y = Y.text, axis.title = a.text, title = t.text) +
  scale_color_manual(values = c("#C3D7A4", "#52854C")) + 
  theme(legend.text=element_text(size=18)) +
  facet_wrap(~ cell_struc,nrow=3) 


########################
# Plot legend
########################
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("topleft", legend =c("Intensity","non-Intensity"), 
       pch=16, pt.cex=3, cex=1.5, bty='n',
       col = c("#C3D7A4", "#52854C"))
mtext("Feature type", at=0.3, cex=2, col = "black")