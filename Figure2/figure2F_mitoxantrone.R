

########################
# R script for figure 2f: Per well density of mitoxantrone treated cells
########################


# clear variables
rm(list=ls())


########################
# load libraries
library(easypackages)
library(ggridges)
libraries("rmarkdown","resample","assertthat","lattice","gplots","ggplot2",
          "fields","readr","data.table","grid","MASS", "reshape", "plyr", "ggpubr", "tictoc")  
libraries("plotly","gridExtra","Rmisc","Hmisc","doBy","caret","corrplot","GGally", "easypackages")  


    ########################
    # import single cell data


    # setwd("~/Directory.../.../")

    setwd("/Users/pearsy02/Documents/NYU_work/Sept_17_2018/R_code_b/Manuscript_Comm_Bio_scripts/Data_files/Fig2")

    cell.dat = read_csv("cell_data_mitoxantrone_figure2f.csv")


    ######################
    title4 = paste("Mitoxantrone", sep = "")
    title4

    title5 = paste("(\u03BC","M)", sep = "")
    title5

    title6 = paste(title4, title5, sep = " ")
    title6



    ######################
    unique(cell.dat$Mitoxantrone_uM)
    count(cell.dat$Mitoxantrone_uM)
    cell.dat$Mitoxantrone_uM = as.character(cell.dat$Mitoxantrone_uM)
    str(cell.dat$Mitoxantrone_uM)

    ######################
    Y.text <- element_text( color = "black", size = 11, angle = 0)
    X.text <- element_text( color = "black", size = 11, angle = 0)
    a.text = element_text(size=15,color = "black")
    t.text = element_text(size=15,color = "black")
    ######################
    
  p.t = 
  ggplot(cell.dat, aes(ObjectTotalIntenCh1, color = Mitoxantrone_uM)) +
  geom_density(size = 0.8) + 
  scale_color_brewer(palette = "Blues", name = title6) +

  theme_classic() +
  theme(legend.position ='none') + 
  theme(axis.text.y = Y.text, axis.text.x = X.text, axis.title = a.text,
        title = t.text, plot.title = element_text(size = 14) ) +
  labs(x = "Total Nucleus Intensity", title = "Dose response", y = NULL) #+

  # adjust axis
  p.t + scale_y_continuous(limits=c(0, .00000020), breaks=c(0, .00000005,
                                                          .00000010, .00000015, .00000020), expand = c(0,0)) +
  scale_x_continuous(limits=c(0, 55297885), breaks=c(0,10000000,20000000, 30000000, 40000000, 50000000,60000000), expand = c(0,0)) 



