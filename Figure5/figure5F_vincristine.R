


########################
# R script for manuscript figure 5 F
# Radial plot of residual scores for vincristine treatment at multiple concentrations, relative to control median
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
# Import 69 feature EMD profile
########################

# setwd("~/Documents/...")

df <- read_csv("ScaledEMDprofile_69feats_785treatments.csv")


########################
# adding metadata here for radial plot figures
########################

df$order = "1"
ind1 = grep("1", df$Metadata_dil_factor)
df$order[ind1] = "1"
ind2 = grep("0.5", df$Metadata_dil_factor)
df$order[ind2] = "2"
ind3 = grep("0.25", df$Metadata_dil_factor)
df$order[ind3] = "3"
ind4 = grep("0.125", df$Metadata_dil_factor)
df$order[ind4] = "4"
ind5 = grep("0.0625", df$Metadata_dil_factor)
df$order[ind5] = "5"
ind6 = grep("0.03125", df$Metadata_dil_factor)
df$order[ind6]="6"
ind7 = grep("0.015625", df$Metadata_dil_factor)
df$order[ind7]="7"
df = df[,c(dim(df)[2],1:dim(df)[2]-1)]

  # Check meta data drug name
  unique(df$Drug_name)
  
  # iidentify vincristine profiles
  drug.current = "vincristine"
  t.keep = which(df$Drug_name %in% drug.current)
  ind4 = grep("control",df$Treatment)
  
  # median of control calculation
  df.c =  df[ind4,]
  c.med = apply(df.c[,22:dim(df.c)[2]], 2, median)
  df.c1 = cbind(df.c[1,1:21],t(c.med))
  df.c1$drug_conc[1] = "DMSO_median"
  
  
  
  ########################
  # order of data - DMSO curves, DMSO median, Treatments
  dat.t = df[t.keep,]
  dat.med = df.c1
  dat.c = df.c
  
  ########################
  # combine all dmso fingerprints with the median of the control
  dat.c1 = rbind(dat.c, dat.med)
  
  ########################
  # annotated data used for later
  df_ch2 = data.frame(rbind(dat.c1, dat.t))
  
  ########################
  # data numeric
  df.num = df_ch2[,21:dim(df_ch2)[2]]
  rownames(df.num) = df_ch2$drug_conc
  ll = length(names(df.num))
  ll
  data <- rbind(rep(1.8,ll) , rep(0,ll) , df.num)
  data.use = data
  
  ###################################
  # dataframe created just for color coding the dose response
  # for annotation and coloring
  
  ind2t = grep(drug.current,df_ch2$Drug_name)
  ind2c = grep("DMSO",df_ch2$Drug_name)
  dat2 = df_ch2[c(ind2c,ind2t),] 
  
  
  ####################################
  # Pick a color palette
  ####################################
  
  # mypalette <- brewer.pal(7,"Blues")
  # image(1:7,1,as.matrix(1:7),col = brewer.pal(7,"Blues"), xlab ="Blues (sequential)",ylab="",xaxt="n",yaxt="n",bty="n")
  
  # mypalette <- brewer.pal(7,"Reds")
  # image(1:7,1,as.matrix(1:7),col=brewer.pal(7,"Reds"), xlab ="Reds (sequential)",ylab="",xaxt="n",yaxt="n",bty="n")
  
  mypalette <- brewer.pal(7,"Purples")
  # image(1:7,1,as.matrix(1:7),col=brewer.pal(7,"Purples"),xlab ="Purples (sequential)",ylab="",xaxt="n",yaxt="n",bty="n")
  
  # mypalette <- brewer.pal(7,"Greens")
  # image(1:7,1,as.matrix(1:7),col = brewer.pal(7,"Greens"), xlab ="Greens (sequential)",ylab="",xaxt="n",yaxt="n",bty="n")
  
  col = mypalette

  
  ####################################
  # adding color metadata 
  ####################################
  
  dat2$cols = col[7]
  dat2$cols[grep("1", dat2$order)]=col[7]
  dat2$cols[grep("2", dat2$order)]=col[6]
  dat2$cols[grep("3", dat2$order)]=col[5]
  dat2$cols[grep("4", dat2$order)]=col[4]
  dat2$cols[grep("5", dat2$order)]=col[3]
  dat2$cols[grep("6", dat2$order)]=col[2]
  dat2$cols[grep("7", dat2$order)]=col[1]
  
  ####################################
  # move last column to first place
  ####################################
  
  dat2 = dat2 [, c(dim(dat2)[2],1:dim(dat2)[2]-1) ]
  
  ord.use = order(dat2$order)
  
  dat2 = dat2[order(dat2$order),]
  
  # call the first two rows from *data.use*
  dat1.a = data.use[1:2,]
  
  # fix this
  dat1.b = data.use[3:dim(data.use)[1],]
  # order it 
  dat1.b = dat1.b[ord.use,]
  # bind the first two rows with this datafrane
  dat1 = rbind(dat1.a, dat1.b)
  
  
  ####################################
  # coloring labels
  ch.col = c("blue", "green", "darksalmon", "deeppink2",
             "tan1", "orange", "red", 
             "turquoise", 
             "orchid", "orchid4")
  ####################################
  
  four_col <- unlist(lapply(names(dat1),function(x){
    
    # panel A
    if (grepl('_NUC_A',x)) ch.col[1] 
    else if (grepl('_RNA',x)) ch.col[2]  
    else if (grepl('_PMG',x)) ch.col[3] 
    else if (grepl('_MITO',x)) ch.col[4]  
    else if (grepl('_NUC_Tex',x)) ch.col[1] 
    
    # panel B
    else if (grepl('_NUC_B',x)) ch.col[1]  
    else if (grepl('_Lyso',x)) ch.col[5]  
    else if (grepl('_Perox',x)) ch.col[6]  
    else if (grepl('_Lipids',x)) ch.col[7] 
    
    # panel C1
    else if (grepl('_NUC_C1',x)) ch.col[1]  
    else if (grepl('_ER',x)) ch.col[8]  
    
    # panel C2
    else if (grepl('_NUC_C2',x)) ch.col[1]  
    else if (grepl('_TUB',x)) ch.col[9] 
    else if (grepl('_ACT',x)) ch.col[10] 
    
    else 'lightgrey'
    
  }))
  
  ###################################
  # customize the radar chart
  ###################################
  
  ###################################
  # remove all the extra DMSO curves
  # ** this is optional **
  
  dat1.sub = dat1[-c(3:315),]
  dat1.sub = dat1[-c(3:332),]
  
  ####################################
  # (1) color of each curve - *pcol*
  colors_border = c("black", dat2$cols[332:338])
  length(colors_border)
  
  # (2) color of interior of fingerprint - *pfcol*
  ins = alpha(rep(NA,8), 0.5)
  length(ins)
  
  # (3) line width specification - *plwd*
  line.width = c(2,rep(1.5,7))
  length(line.width)
  
  # (4) line type specification - *plty*
  line.type = c(2,rep(1,7))
  length(line.type)
  
  # (5) line shapes specification - *pty*
  line.shape = c(rep(NA,8))
  length(line.shape)
  
  # check dimensions again
  dim(dat1)
  

  
  ###################################
  # The radial plot
  ###################################
  
  title.use = paste("Compound", drug.current, sep = ": ")
  
  radarchart2(dat1.sub, axistype = 4 , seg = 2,
              
              #custom polygon
              pcol = colors_border , pfcol = ins, plwd = line.width , plty = line.type, pty = line.shape,
              #custom the grid
              cglcol="grey", cglty = 1, axislabcol="grey",  cglwd=.6, 
              #custom labels
              vlcex = 0.5, calcex = 1, palcex = 1, vlabcol = four_col,
              title = title.use
  )
  
   
  
  ###################################
  # Expand the range in the data for better visualization
  ###################################
  
  dat2.sub = dat1.sub + 0.5
  dat2.sub[1,] = 1.5
  dat2.sub[2,] = 0
  
  
  ##################
  # Axis seg = 3
  ##################
  
  radarchart2(dat2.sub, axistype = 4 , seg = 3,
              
              #custom polygon
              pcol = colors_border , pfcol = ins, plwd = line.width , plty = line.type, pty = line.shape,
              #custom the grid
              cglcol="grey", cglty=1, axislabcol="grey",   cglwd=.6,
              #custom labels
              vlcex = 0.5, calcex = 1, palcex = 1, vlabcol = four_col,
              title = title.use
  )
  
  
  
  
  ###################################
  # Residual radial plots - standardized to the control
  ##################################
  
  num.f = dim(dat1.sub)[2]
  dat.resid = dat1.sub
  
  for (i in seq(1,num.f) ) {
  
    temp = dat1.sub[,i] - dat1.sub[3,i]
  
    dat.resid[,i] = temp
  
}
    
  dat.resid[1,] = 1
  dat.resid[2,] = 0
  
  
  radarchart2(dat.resid, axistype = 4 , seg = 2,
              
              #custom polygon
              pcol = colors_border , pfcol = ins, plwd = line.width , plty = line.type, pty = line.shape,
              #custom the grid
              cglcol="grey", cglty=1, axislabcol="grey",   cglwd=.6,
              #custom labels
              vlcex = 0.5, calcex = 1, palcex = 1, vlabcol = four_col,
              title = title.use
  )
  
  
  ###################################
  # Residual radial plots - standardized to the control
  # Expand the range for visualization
  ##################################
  
  dat.resid = dat.resid + 0.5
  dat.resid[1,] = 1
  dat.resid[2,] = 0
  
  radarchart2(dat.resid, axistype = 4 , seg = 2,
              
              #custom polygon
              pcol = colors_border , pfcol = ins, plwd = line.width , plty = line.type, pty = line.shape,
              #custom the grid
              cglcol="grey", cglty=1, axislabcol="grey",   cglwd=.6,
              #custom labels
              vlcex = 0.5, calcex = 1, palcex = 1, vlabcol = four_col,
              title = title.use
  )
  
  
  
  
  ############
  # legend for vincristine - start
  ############ 
  
  title4 = paste("Vincristine", sep = "")
  title4
  
  title5 = paste("(\u03BC","M)", sep = "")
  title5
  
  title6 = paste(title4, title5, sep = " ")
  title6

  plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
  
  mypalette<-brewer.pal(7,"Purples")
  # image(1:7,1,as.matrix(1:7),col=mypalette,xlab="Reds (sequential)",ylab="",xaxt="n",yaxt="n",bty="n")
  col.tox = mypalette
  
  
  unique(df1$Metadata_dil_factor)
  # [1] "0.25"     "0.0625"   "0.03125"  "0.015625" "0.5"      "1"        "0.125"    "NULL"  
  labels = c("0.3125", "0.625", "1.25", "2.5", "5.0","10.0", "20.0")
  
  legend("center", 
         legend = labels, 
         col = col.tox,
         #pch = c(19,19, 19, 19, 19, 19, 19, 19, 19,19 ), 
         pch = c(15, 15, 15, 15, 15, 15, 15, 15, 15, 15 ), 
         # pch = c(16, 16, 16, 16, 16, 16, 16, 16, 16, 16 ), 
         
         bty = "n", 
         pt.cex = 2, 
         cex = 1.2, 
         text.col = "black", 
         horiz = F , 
         inset = c(0.1, 0.1))
  
  mtext(title6, at=1.0, cex=1.5)
  
  ############
  # legends - end
  ############ 
  
  
  
  
  