


########################
# R script for figure 2E: Cell counts heatmap
########################


########################
# load libraries
########################

library(easypackages)
library(resample) 
libraries("assertthat","lattice","gplots","ggplot2","fields","readr","data.table","grid","MASS")  
libraries("plotly","gridExtra","Rmisc","Hmisc","doBy","caret","corrplot","GGally", "easypackages")  
library(reshape) 
library(plyr)
library(ggpubr)
library(RColorBrewer)

########################
# import the annotated well median data
########################


well_med_dat <- # raw_medians_fig2E.csv
  
names(well_med_dat)
dat =  well_med_dat
temp1 = dat
 
######################
# identify outlier wells and empty wells

empty = which(grepl("empty", dat$Treatment)) 
low.c = which(dat$CELL_COUNT < 40)

# set those to NA
temp1[empty,13:dim(temp1)[2]] = NA
temp1[low.c,13:dim(temp1)[2]] = NA


######################
# saving heatmap matrices in "lists"

datalist.plates = list()
datalist.plates.all = list()
datalist.plates.c = list()
datalist.plates.all.c = list()
merged.save=list()


panelA.inten = temp1
panelA.inten.num = panelA.inten[,13:dim(panelA.inten)[2]]
plate.unique = sort(unique(panelA.inten$PlateNumber))

  # j is a feature
  j = 1
  # loop through each plate
  for (i in seq(1,length(plate.unique))) {     
    
    plate.ind = which(panelA.inten$PlateNumber == plate.unique[i])
    plate.current = panelA.inten.num[plate.ind,j]
    plate.current.labels = panelA.inten[plate.ind,]
    plate.current.ord <- plate.current.labels[order(plate.current.labels$WellId),]
    
    ######################
    # The data
    plate.use <- plate.current.ord[,12+j]
    
    ######################
    # The feature name
    feat.use <- names(plate.current.ord)[12+j]
    feat.name <- names(panelA.inten.num)[j]
    
    ######################
    # raw data matrix
    feature1_mat <- data.matrix(plate.use)
    feature1_mat1 <- matrix(feature1_mat, nrow = 14, ncol = 22, byrow = TRUE)
    rownames(feature1_mat1) <- LETTERS[2:15]
    colnames(feature1_mat1) <- 2:23
    
    ######################
    # full matrix save
    datalist.plates[[i]] <- feature1_mat1 
    
    ######################
    # heatmap color
    coul <- colorRampPalette(brewer.pal(8, "RdBu"))(25)
    
    ######################
    # empty treatment wells set to NA
    ind.t = which(grepl("treated",plate.current.ord$Treatment)|grepl("empty",plate.current.ord$Treatment))
    ind.c = which(grepl("control",plate.current.ord$Treatment)|grepl("empty",plate.current.ord$Treatment))

    
    plate.use.c = plate.current.ord[,12+j] # setting all treatment wells to NA
    plate.use.c[ind.t,] = NA
    feature1_mat.ct <- data.matrix(plate.use.c)
    feature1_mat.ct1 <- matrix(feature1_mat.ct, nrow = 14, ncol = 22, byrow = TRUE) # raw data
    rownames(feature1_mat.ct1) <- LETTERS[2:15]
    colnames(feature1_mat.ct1) <- 2:23
    
    # save full matrix - control
    datalist.plates.c[[i]] <- feature1_mat.ct1 
    
  }
  
  
  
  ######################
  # binding and storing all matricies
  
  A = rbind(datalist.plates[[1]],datalist.plates[[2]], datalist.plates[[3]])
  B = rbind(datalist.plates[[4]],datalist.plates[[5]], datalist.plates[[6]])
  C = cbind(A,B)
  datalist.plates.all[[j]] <- C
  full = C
  
  A = rbind(datalist.plates.c[[1]],datalist.plates.c[[2]], datalist.plates.c[[3]])
  B = rbind(datalist.plates.c[[4]],datalist.plates.c[[5]], datalist.plates.c[[6]])
  C = cbind(A,B)
  datalist.plates.all.c[[j]] <- C
  controls = C

  merged = cbind(full,controls)
  merged.save[[j]] = merged
  
  ######################
  # features
  feat.name <- names(panelA.inten.num) 

  ######################
  # heatmap colors 
  new.palette <- colorRampPalette(brewer.pal(8, "RdBu"))(25)

  ######################
  # count features
    i = 1

  ######################
  # the heatmap plot
  
  levelplot(merged.save[[i]],panel = function(...){
    panel.levelplot(...)
    
    panel.abline(h = 22.5, col = "black", lwd = .5)
    panel.abline(h = 44.5, col = "black", lwd = .5)
    panel.abline(h = 67.5, col = "black", lwd = .5)
    
    panel.abline(v = 14.5, col = "black", lwd = .5)
    panel.abline(v = 28.5, col = "black", lwd = .5)
  },
  # at = seq(rg1[1], rg1[2], inc1),
  col.regions= new.palette,
  xlab =  paste(feat.name[i]),
  ylab = NULL,
  main = NULL,
  scale=list(x=list(at=c(1,14), rot=90, cex=0.8), y=list(at=c(1,22), rot=90, cex=0.8)), aspect = "fill")
  

  


