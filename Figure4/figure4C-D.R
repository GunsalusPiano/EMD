

########################
# R script for figure  4 C & D: Distribution of statistical distances
########################



########################
# Load libraries
########################
library(easypackages);library(resample) ;#library(caroline)
libraries("assertthat","lattice","gplots","ggplot2","fields","readr","data.table","grid","MASS")  
libraries("plotly","gridExtra","Rmisc","Hmisc","doBy","caret","corrplot","GGally", "easypackages")  
library(reshape);
library(plyr);library(ggpubr);
library("scatterplot3d") ;library(rgl);library(rafalib);library(ggridges)#;library("FactoClass")


########################
# Import output from the replicate reproducibility analysis   
########################

reps.dmso <- # control_data_fig4C.csv 
reps.trt <- # treatment_data_fig4C.csv

########################
# check data variables
unique(reps.dmso$variable)

reps.dmso$treatment = "DMSO"
reps.trt$treatment = "TRT"

########################
# melt the data 
df.melt = rbind(reps.dmso, reps.trt)
df.melt$variable_treatment = paste(df.melt$variable, df.melt$treatment, sep = "_")


######################
# ggplot parameters
######################

Y.text <- element_text( color = "black", size = 15, angle = 0)
X.text <- element_text( color = "black",  size = 15, angle = 0)
a.text = element_text(size=20, color = "black") 
t.text = element_text(size=17, color = "black")
strip.text.a = element_text(size = 8, color = "black")

######################
# Manuscript figure 4 C
######################

ggplot(df.melt, aes(x = value, color = variable)) +
  geom_density(size = .8, aes(linetype = treatment)) + 
  theme_bw() +
  scale_color_manual(name = "Statistic", values=c('ks_dist' = "#F21A00", 'clip_emd_dist' = "purple",  'clip_z_dist' = "cyan")) +
  scale_linetype_manual(name = "Treatment", values=c('DMSO' = "solid",'TRT' = "dashed")) +
  theme(axis.text.y = Y.text, axis.text.x = X.text, axis.title = a.text, title = t.text, plot.title = element_text(size=18) ) +
  theme(legend.position='right') + 
  labs(x = NULL, y = "Density", title = NULL) # +
# facet_wrap(~treatment)







######################
# Analysis for manuscript figure 4 D
# Ranking features based on reproducibility performance 
######################



######################
# Import data files from replicate reproducibility testing
# Includes metadata such as plate_well_A and plate_well_B - these columns show which wells were compared
######################

# DMSO replicates 
reps.dmso1 <- # control_replicate_test_allpanels_fig4D.csv

# Treatment replicates
reps.trt1 <- # treatment_replicate_test_allpanels_fig4D.csv

  


######################
# aggregate - find mean of each feature - DMSO replicates
######################

meds.ks = aggregate(reps.dmso1[, 10], list(reps.dmso1$feature), mean, na.rm = TRUE)
meds.emd = aggregate(reps.dmso1[, 11], list(reps.dmso1$feature), mean, na.rm = TRUE)

meds = cbind(meds.ks, meds.emd)
meds.df = meds[,c(1,2,4)]
names(meds.df)[1] = "feature"
meds.dmso = meds.df
meds.dmso$emd_dist_dmso = meds.dmso$emd_dist

######################
# aggregate - find mean of each feature - treatment replicates
######################

meds.ks = aggregate(reps.trt1[, 10], list(reps.trt1$feature), mean, na.rm = TRUE)
meds.emd = aggregate(reps.trt1[, 11], list(reps.trt1$feature), mean, na.rm = TRUE)

meds = cbind(meds.ks, meds.emd)

meds.df = meds[,c(1,2,4)]
names(meds.df)[1] = "feature"
meds.trt = meds.df
meds.trt$emd_dist_trt = meds.trt$emd_dist


######################
# combine the treatment and DMSO-control per feature medians
######################

emd.all = cbind(meds.dmso, meds.trt)

######################
# EMD first
data = emd.all$emd_dist_dmso
######################

hist(data)
plot(data)
median(data)

lowerq = quantile(data)[2]
upperq = quantile(data)[4]
iqr = upperq - lowerq 
iqr

######################
# 75%  
mild.threshold.upper = (iqr * 1.5) + upperq    
mild.threshold.upper # 0.6521563
mild.threshold.lower = lowerq - (iqr * 1.5)
######################

df = emd.all
df1=df[order(df$emd_dist_dmso),]
df1$ind=1:174
df1$ind_feature=paste(df1$ind,df1$feature, sep = "_")
df1$feat_ord = as.factor(df1$feature)
df1$feat_ord <- factor(df1$feat_ord, levels = df1$feat_ord[order(df1$ind)])
df1$feat_ord  # notice the changed order of factor levels

######################
# EMD sorted by mean plot TRT and DMSO
######################

df1 = df1[,c(1:4,8:11)]
df1$noise = "keep"
df1$noise[160:174] = "drop"

a <- ifelse(df1$noise == "drop", "red", "black")


######################
# Plot parameters
X.text <- element_text(family = "sans", size =7, angle = 65, hjust=1, color = a)
Y.text <- element_text(color = "black", size = 14, angle = 0)

######################
# Plot of sorted average EMD values for each feature

p = ggplot(df1, aes(x = feat_ord, y = emd_dist_dmso)) +

    theme_minimal() +
    geom_point(size = 1) +
    geom_point(data = df1, mapping=aes(x = feat_ord, y = emd_dist_trt), colour = "gray", size= 1.5, shape = 17) +
    geom_vline(aes(xintercept=159), color = "red", linetype = "dashed", size = .6) +
    geom_hline(aes(yintercept = mild.threshold.upper), color = "red", linetype = "dashed", size = .6) +
    theme(strip.text = strip.text.a, axis.text.y = Y.text, axis.text.x = X.text, 
        axis.title = a.text, title = t.text, plot.title = element_text(size=22) ) +
    labs(y = "Average EMD score",  x = NULL , title = NULL)

######################
# Manually adding Y axis scale
######################
p + scale_y_continuous(breaks=seq(0,3,.5))