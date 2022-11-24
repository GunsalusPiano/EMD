
########################
# R script for Supplementary figure  4  and manuscript figure 5
# Figure 5. D & E: Replicate aggregation and EMD profiling using global controls.
# Supplementary Figure 4: Full feature EMD fingerprints of individual control samples.
########################


########################
# Clear variables
########################

rm(list=ls())

########################
# Load libraries
########################

library(car)
library(easypackages)
libraries("assertthat","lattice","gplots","ggplot2","fields","readr","data.table","grid","MASS")  
libraries("plotly","gridExtra","Rmisc","Hmisc","doBy","caret","corrplot","GGally", "easypackages")  
library(viridis)
library(RColorBrewer)
library(ggplot2)
library(ggrepel)
library("FactoMineR")
library("factoextra")
library(fmsb)

########################
# Run the function (radarchart2) for plotting below
########################

########################
# set working directory

# setwd("~/Documents/...")
source("radarchart2_new.R")


########################
# Import full 174 feature profile
########################

# setwd("~/Documents/...")
df <- read_csv("raw_emd_profile174.csv")

########################
# Keep 5 main meta data columns
########################

st = grep("ObjectArea_NUC_A", names(df))
df_ch1 = df[,c(1,8,9,6,10,11, st:dim(df)[2])]

########################
# Scale the full feature profile
########################

########################
# Log transform the feature data again
########################

dlog = log(df_ch1[,7:180]) # natural log

max.log <- max(dlog, na.rm = TRUE)
max.log

min.log <- min(dlog, na.rm = TRUE)
min.log

dlog.sc = (dlog-min.log)/(max.log - min.log)

########################
# combine numeric dataframe with metadata
d.sc = cbind(df_ch1[,1:6], dlog.sc)
d.sc.temp = d.sc
df_ch1 = d.sc

########################
# Adding metadata
########################

df_ch1$order = "1"
ind1=grep("1", df_ch1$Metadata_dil_factor)
df_ch1$order[ind1]="1"
ind2=grep("0.5", df_ch1$Metadata_dil_factor)
df_ch1$order[ind2]="2"
ind3=grep("0.25", df_ch1$Metadata_dil_factor)
df_ch1$order[ind3]="3"
ind4=grep("0.125", df_ch1$Metadata_dil_factor)
df_ch1$order[ind4]="4"
ind5=grep("0.0625", df_ch1$Metadata_dil_factor)
df_ch1$order[ind5]="5"
ind6=grep("0.03125", df_ch1$Metadata_dil_factor)
df_ch1$order[ind6]="6"
ind7=grep("0.015625", df_ch1$Metadata_dil_factor)
df_ch1$order[ind7]="7"

df_ch1 = df_ch1[,c(dim(df_ch1)[2],1:dim(df_ch1)[2]-1)]

########################
# Identify control samples
########################

t.keep = which(df_ch1$Metadata_Drug_name %in% "DMSO")
ind4 = grep("control", df_ch1$Metadata_type)

########################
# Median of control samples
########################

df.c =  df_ch1[ind4,]
c.med = apply(df.c[,8:dim(df.c)[2]], 2, median)
df.c1 = cbind(df.c[1,1:7],t(c.med))
df.c1$drug_conc = "DMSO_median"

########################
# Data preparation for radial plot
########################

df_ch2 = data.frame(rbind(df_ch1[t.keep,], df.c1))
df.num = df_ch2[,8:dim(df_ch2)[2]]
rownames(df.num) = df_ch2$drug_conc

ll = length(names(df.num))
data <- rbind(rep(1,ll) , rep(0,ll) , df.num)
data.use = data

########################
# Rows that need to be preserved for the radial plot
########################

rownames(data.use)[c(1,2,length(rownames(data.use)))]
dat.rad = data.use[c(1,2,length(rownames(data.use))),]
kp = c(1,2,length(rownames(data.use)))

ind1 = grep("_DMSO",rownames(data.use))
dat1 = data.use[c(kp, ind1),]


########################
# customize the radatchar
########################

dat1 = rbind(dat1,dat1[3,])

########################
# Coloring of DMSO fingerprints
########################

colors_border = c("black",rep("lightgrey", 330), "black")

ins = alpha(c(NA,rep(NA,330)), 0.5)

########################
# Coloring of labels
########################

col.one1 = c("blue", "green", "darksalmon", "deeppink2",
             "tan1", "orange", "red", 
             "turquoise", 
             "orchid", "orchid4")

channel_colors1 <- unlist(lapply(names(dat1),function(x){
  
  if(grepl("NUC_A",x))  col.one1[1]        
  else if (grepl("RNA",x)) col.one1[2]     
  else if (grepl("PMG",x))  col.one1[3]    
  else if (grepl("MITO",x)) col.one1[4]    
  else if (grepl("NUC_Tex",x)) col.one1[1] 
  
  else if (grepl("NUC_B",x)) col.one1[1]   
  else if (grepl("Lyso",x)) col.one1[5]    
  else if (grepl("Perox",x)) col.one1[6]   
  else if (grepl("Lipids",x)) col.one1[7]  
  
  else if (grepl("NUC_C1",x)) col.one1[1]  
  else if (grepl("ER",x)) col.one1[8]     
  
  else if (grepl("NUC_C2",x)) col.one1[1]  
  else if (grepl("TUB",x))    col.one1[9] 
  else if (grepl("ACT",x))    col.one1[10] 
  else 'lightgrey'
  
  
}))





########################
# ‘pcol’ → line color
# ‘pfcol’ → fill color
# ‘plwd’ → line width
########################

#########################
# *** Radial plot of MEDIAN of controls ***
#########################

dat.f = dat1

radarchart2(dat.f[1:3,], axistype = 4 , seg = 2,
            
            #custom polygon
            pcol = colors_border , pfcol = ins, plwd = c(1.5,rep(.5,330),1.5) , plty = c(5,rep(1,330),5), pty = c(rep(NA,330)),
            
            # caxislabels = seq(from = 0, to = 1, length = 4),
            
            #custom the grid
            cglcol = "lightgrey", cglty = 1, axislabcol = "lightgrey",    cglwd = .4,
            
            # Custom labels # ‘vlcex’ → group labels size
            vlcex = 0.5, calcex = 1, palcex = 1, vlabcol = channel_colors1,
            title = paste("DMSO samples")
)


##############################
# All control fingerprints
# axistype = 6 - no labels
# axistype = 4 - labels
##############################

# *** CAUTION: note that with individual control fingerprints the plot takes a long time ***

radarchart2(dat.f, axistype = 4 , seg = 2,
            
            #custom polygon
            pcol = colors_border , pfcol = ins, plwd = c(1.5,rep(.5,330),1.5) , plty = c(5,rep(1,330),5), pty = c(rep(NA,330)),
            
            # caxislabels = seq(from = 0, to = 1, length = 4),
            
            #custom the grid
            cglcol = "lightgrey", cglty = 1, axislabcol = "lightgrey",    cglwd = .4,
            
            # Custom labels # ‘vlcex’ → group labels size
            vlcex = 0.5, calcex = 1, palcex = 1, vlabcol = channel_colors1,
            title = paste("DMSO samples")
)







########################
# Next step is to standardize each EMD to the control EMD
########################

########################
# separate treatments from control
ind4 = grep("control", df_ch1$Metadata_type)
df.c =  df_ch1[ind4,]
df.t =  df_ch1[-ind4,]

########################
# median calculation
c.med = apply(df.c[,8:dim(df.c)[2]], 2, median)
c.mad = apply(df.c[,8:dim(df.c)[2]], 2, mad)
c.mean = apply(df.c[,8:dim(df.c)[2]], 2, mean)


########################
#control
dat.num.c = df.c[,8:dim(df.c)[2]]

#treatment
dat.num.t = df.t[,8:dim(df.t)[2]]

#all data
dat.num.all = df_ch1[,8:dim(df_ch1)[2]]

########################
# subtract median of the contol from all values
########################

z.med = sweep(dat.num.all, 2, c.med, '-')

########################
# add the annotations back
emd.residual <- cbind(df_ch1[,1:7],z.med) 

########################
# Check values of the control
########################

ind.c = which(emd.residual$Metadata_type == "control")
min(z.med[ind.c,])  
max(z.med[ind.c,])  

residual.c = z.med[ind.c,]
residual.med = apply(residual.c, 2, median)


########################
# Identify control samples
########################

df_ch1 = emd.residual
t.keep = which(df_ch1$Metadata_Drug_name %in% "DMSO")
ind4 = grep("control", df_ch1$Metadata_type)

########################
# Median of control samples
########################

df.c =  df_ch1[ind4,]
c.med = apply(df.c[,8:dim(df.c)[2]], 2, median)
df.c1 = cbind(df.c[1,1:7],t(c.med))
df.c1$drug_conc = "DMSO_median"

########################
# Data preparation for radial plot
########################

df_ch2 = data.frame(rbind(df_ch1[t.keep,], df.c1))
df.num = df_ch2[,8:dim(df_ch2)[2]]
rownames(df.num) = df_ch2$drug_conc

ll = length(names(df.num))
data <- rbind(rep(1,ll) , rep(0,ll) , df.num)
data.use = data

########################
# Rows that need to be preserved for the radial plot
########################

rownames(data.use)[c(1,2,length(rownames(data.use)))]
dat.rad = data.use[c(1,2,length(rownames(data.use))),]
kp = c(1,2,length(rownames(data.use)))

ind1 = grep("_DMSO",rownames(data.use))
dat1 = data.use[c(kp, ind1),]



########################
# Radial plot of residuals
########################

dat1 = rbind(dat1,dat1[3,])

########################
# Coloring of DMSO fingerprints
########################

colors_border = c("black",rep("lightgrey", 330), "black")

ins = alpha(c(NA,rep(NA,330)), 0.5)

########################
# Coloring of labels
########################

col.one1 = c("blue", "green", "darksalmon", "deeppink2",
             "tan1", "orange", "red", 
             "turquoise", 
             "orchid", "orchid4")

channel_colors1 <- unlist(lapply(names(dat1),function(x){
  
  if(grepl("NUC_A",x))  col.one1[1]        
  else if (grepl("RNA",x)) col.one1[2]     
  else if (grepl("PMG",x))  col.one1[3]    
  else if (grepl("MITO",x)) col.one1[4]    
  else if (grepl("NUC_Tex",x)) col.one1[1] 
  
  else if (grepl("NUC_B",x)) col.one1[1]   
  else if (grepl("Lyso",x)) col.one1[5]    
  else if (grepl("Perox",x)) col.one1[6]   
  else if (grepl("Lipids",x)) col.one1[7]  
  
  else if (grepl("NUC_C1",x)) col.one1[1]  
  else if (grepl("ER",x)) col.one1[8]     
  
  else if (grepl("NUC_C2",x)) col.one1[1]  
  else if (grepl("TUB",x))    col.one1[9] 
  else if (grepl("ACT",x))    col.one1[10] 
  else 'lightgrey'
  
  
}))

##############################
# Median of the control fingerprint
# axistype = 6 - no labels
# axistype = 4 - labels
##############################

radarchart2(dat1[1:3,], axistype = 4 , seg = 2,
            
            #custom polygon
            pcol = colors_border , pfcol = ins, plwd=c(2,rep(.5,313),2) , plty = c(2,rep(1,313),2), pty = c(rep(NA,315)),
            #custom the grid
            cglcol="grey", cglty=1, axislabcol="grey",    cglwd=.4,
            #custom labels
            vlcex = 0.5, calcex = 1, palcex = 1, vlabcol = channel_colors1,
            title = paste("DMSO samples")
)

##############################
# Scale up the values for better visualization
##############################

dat.positive = dat1 + 0.5
dat.positive[1,] = 1
dat.positive[2,] = 0

# *** CAUTION: note that with individual control fingerprints the plot takes a long time ***

# Number of control fingerprints to plot
N = 10

radarchart2(dat.positive[c(1:N,334),], axistype = 4 , seg = 2,
            
            #custom polygon
            pcol = colors_border , pfcol = ins, plwd=c(1.5,rep(.5,330),1.5) , plty = c(5,rep(1,330),5), pty = c(rep(NA,330)),
            #custom the grid
            cglcol="lightgrey", cglty=1, axislabcol="lightgrey",    cglwd=.4,
            #custom labels
            vlcex = 0.5, calcex = 1, palcex = 1, vlabcol = channel_colors1,
            title = paste("Residuals")
)







