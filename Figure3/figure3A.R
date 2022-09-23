

  ########################
  # R script for figure 3A: Detection of non-uniformity among control wells
  ########################



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


  ########################
  # Import well level (median) panel A data from six plates   
  ########################
  
  # well_medians_rawdata_panelA.csv

   
  ########################
  # Remove non-biological features - mito circ and nucleus ring features
  ########################
    
  dat =  panelA.med
  rr = names(dat)
  ind <- which(grepl("Ring", rr) & grepl("NUC_Tex", rr)|grepl("Circ", rr) & grepl("MITO", rr))
  rr[ind]
  dat1 = dat[,-(ind)] 

  ########################
  # Remove empty or low-count wells
  ########################

  empty = which(grepl("empty", dat1$Treatment)) 
  low.c = which(dat1$CELL_COUNT < 40)

  ########################
  # Save full data set to variable "dat2"
  ########################
  
  dat2 = dat1

  ########################
  # Set feature data to NA for the case of noisy wells/low count/no cells
  ########################

  dat1[empty,13:dim(dat1)[2]] = NA
  dat1[low.c,13:dim(dat1)[2]] = NA


  ########################
  # Lists to use for saving matrices
  ########################

  ########################
  datalist.plates = list()
  datalist.plates.all = list()
  ########################
  datalist.plates.c = list()
  datalist.plates.all.c = list()
  ########################
  
  merged.save=list()



  ########################
  # For each feature and each plate
  ########################

  panelA.inten = dat1
  panelA.inten.num = panelA.inten[,13:dim(panelA.inten)[2]]

  ########################
  # Identify each plate
  ########################
  
  plate.unique = sort(unique(panelA.inten$PlateNumber))

  ########################
  # Parse through each feature (j) and each plate (j)
  ########################
  
  for (j in seq(1,length(panelA.inten.num))) { 
  
  for (i in seq(1,length(plate.unique))) {   
    
    plate.ind = which(panelA.inten$PlateNumber == plate.unique[i])
    plate.current = panelA.inten.num[plate.ind,j]
    plate.current.labels = panelA.inten[plate.ind,]
    plate.current.ord <- plate.current.labels[order(plate.current.labels$WellId),]
    
    ########################
    # The data/feature
    ########################
    plate.use <- plate.current.ord[,12+j]
    
    ########################
    # The feature name
    ########################
    feat.use <- names(plate.current.ord)[12+j]
    feat.name <- names(panelA.inten.num)[j]
    
    ########################
    # Raw data matrix
    ########################
    
    feature1_mat <- data.matrix(plate.use)
    feature1_mat1 <- matrix(feature1_mat, nrow = 14, ncol = 22, byrow = TRUE)
    rownames(feature1_mat1) <- LETTERS[2:15]
    colnames(feature1_mat1) <- 2:23
    
    ########################
    # Saving raw data matrix to list
    ########################
    
    datalist.plates[[i]] <- feature1_mat1 

    
    ########################
    # Indices of treatment wells
    # Setting treatments/empty wells to NA
    ########################
    
    ind.t = which(grepl("treated",plate.current.ord$Treatment)|grepl("empty",plate.current.ord$Treatment))
    ind.c = which(grepl("control",plate.current.ord$Treatment)|grepl("empty",plate.current.ord$Treatment))
    
    plate.use.c = plate.current.ord[,12+j] 
    
    ########################
    # setting all treatment wells to NA
    plate.use.c[ind.t,] = NA
    
    ########################
    # formatting the matrix/heatmap to copy the plate layout (rows, columns and well positions)
    feature1_mat.ct <- data.matrix(plate.use.c)
    feature1_mat.ct1 <- matrix(feature1_mat.ct, nrow = 14, ncol = 22, byrow = TRUE) 
    rownames(feature1_mat.ct1) <- LETTERS[2:15]
    colnames(feature1_mat.ct1) <- 2:23
    
    
    ########################
    # saving the dmso control matrix to a list
    datalist.plates.c[[i]] <- feature1_mat.ct1 
    
  }
  
  ########################
  # For each feature - all six data matrices/plates data are stored below into new formatted lists
  ########################
  
  A = rbind(datalist.plates[[1]],datalist.plates[[2]], datalist.plates[[3]])
  B = rbind(datalist.plates[[4]],datalist.plates[[5]], datalist.plates[[6]])
  C = cbind(A,B)
  datalist.plates.all[[j]] <- C
  full=C
  
  ########################
  # Formatting the plates, but showing only controls (diagonal wells)
  ########################
  
  A = rbind(datalist.plates.c[[1]],datalist.plates.c[[2]], datalist.plates.c[[3]])
  B = rbind(datalist.plates.c[[4]],datalist.plates.c[[5]], datalist.plates.c[[6]])
  C = cbind(A,B)
  datalist.plates.all.c[[j]] <- C
  controls = C

  ########################
  # Merging plates next to eachother, forming one larger matrix
  
  merged=cbind(full,controls)
  merged.save[[j]] = merged
  
}



  ########################
  # Plotting heatmaps
  ########################

  ########################
  # Identify feature
  feat.name <- names(panelA.inten.num) 

  ########################
  # Define color palette
  new.palette = colorRampPalette(c("black","red","yellow","white"), space="rgb")

  ########################
  # Other color options
  cols <- brewer.pal(3, "BuGn")
  pal <- colorRampPalette(cols)

  ########################
  # Initiate list for saving
  f1.c = list()

  
  # i = 5 -- "ObjectTotalInten_NUC_A"
  
  for (i in  seq(1,length(feat.name))) {  

  ########################
  # Raw data heatmap - all size plates
  
  full <-levelplot(datalist.plates.all[[i]],panel = function(...){
    panel.levelplot(...)
    panel.abline(h = 22.5, col = "black", lwd = .5)
    panel.abline(v = 14.5, col = "black", lwd = .5)
    panel.abline(v = 28.5, col = "black", lwd = .5)

  },
  #at = seq(rg1[1], rg1[2], inc1),
  col.regions= new.palette,
  xlab = paste(feat.name[i]),
  ylab = "RAW DATA",
  main = NULL,
  scale=list(x=list(at=c(1,14), rot=90, cex=0.7), y=list(at=c(1,22), rot=90, cex=0.7)), aspect = "fill")
  
  ########################
  # print heatmap of all wells to R plots
  full
  
  ########################
  # Raw data heatmap of only controls - all plates
  
  controls <-levelplot(datalist.plates.all.c[[i]],panel = function(...){
    panel.levelplot(...)
    panel.abline(h = 22.5, col = "black", lwd = .5)
    panel.abline(v = 14.5, col = "black", lwd = .5)
    panel.abline(v = 28.5, col = "black", lwd = .5)
    
  },
  # at = seq(rg1[1], rg1[2], inc1),
  col.regions= new.palette,
  xlab = paste(feat.name[i]),
  ylab = "Control wells",
  main = NULL,
  scale=list(x=list(at=c(1,14), rot=90, cex=0.7), y=list(at=c(1,22), rot=90, cex=0.7)), aspect = "fill")
  
  ########################
  # print heatmap of control wells to R plots
  controls
  
  ########################
  # Full plate of data on the left, control wells only, on the right
  
  f1.c[[i]] <-levelplot(merged.save[[i]],panel = function(...){
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

  ########################
  # print full plate and control heatmaps side by side to R plots
  f1.c

  # grid.arrange(f1.t, f1.a, f1.bscore2, f1, f1.mp, ncol = 3)
  
}






