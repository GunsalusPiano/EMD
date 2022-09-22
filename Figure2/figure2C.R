

######################
# load libraries
library(easypackages)
libraries("assertthat","lattice","gplots","ggplot2","fields","readr","data.table","grid","MASS")  
libraries("plotly","gridExtra","Rmisc","Hmisc","doBy","caret","corrplot","GGally", "easypackages")  


######################
# read in data matrix
# my_data = # read in csv file

temp1 = my_data
rownames(temp1) = my_data$X1
temp2 = temp1[,-1]
rownames(temp2) = my_data$X1


######################
# convert temp2 to a matrix

temp2_numeric = apply(temp2, 2, as.numeric)
rownames(temp2_numeric)=rownames(temp2)
temp3 = temp2_numeric

######################
# create heatmap of dataframe called "temp3"
heatmap.2(temp3, Colv = T, Rowv = T, trace = "none", density ="none", col = bluered(100), cexRow=.5, cexCol=0.5, 
          
          scale="none",  hclust=function(x) hclust(x,method="average"))

