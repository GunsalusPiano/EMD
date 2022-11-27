########################
# R script for figure  7 C: Cell cycle distributions  
########################


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
library(tsne)
library("scatterplot3d") 
library(rgl)
library(rafalib)
library(tictoc)
library(ggridges)
library(plyr)
library(RColorBrewer)
library(umap)



# set working directory

# setwd("~/Directory.../.../")

# treatments
dat.sub <- read_csv("figure7c_treatments.csv")

# control
dat.c <- read_csv("figure7c_control.csv")


# Four compounds
grp1 = c("tolfenamic-acid")    # green group
grp2 = c("methotrexate") # purple group
grp3 = c("irinotecan") # green group
grp4 = c("tanespimycin") # blue group

######
# grp1
######

in.y <- which(dat.sub$Drug_name %in% grp1)
dat.sub1 = dat.sub[in.y,]

dat.sub1$Drug_name1 = grp1
dat.c$Drug_name1 = grp1

##########
# Plot parameters
##########

a.text = element_text(size=25, color = "blue")
t.text = element_text(size=25, color = "blue")
Y.text <- element_text( color = "black", size = 19, angle = 0)
X.text <- element_text( color = "black", size = 19, angle = 0)
# strip.text.a = element_text(size = 10, color = "blue")
strip.text.a = element_text(size = 24, color = "black")

dat.sub1$Metadata_dil_factor1 = as.character(dat.sub1$Metadata_dil_factor)


cc = ggplot(dat.sub1, aes(ObjectTotalInten_NUC_A, colour = Metadata_dil_factor1)) +
  geom_density(aes(linetype = Treatment, size = Treatment)) +
  geom_density(data = dat.c, aes(ObjectTotalInten_NUC_A, size = Treatment), colour = "black", linetype = "dashed") +
  
  ##########################
  # mapped variables - linetype, size and color
  scale_linetype_manual(values=c("solid", "solid")) +
  scale_size_manual( values = c(.8,.9) ) +
  ##########################
  #colors
  scale_color_brewer(palette = "Greens") + 
  ##########################
  theme_minimal() +
  theme(legend.position = 'none') +
  theme(strip.text = strip.text.a, axis.text.y = Y.text, 
        axis.text.x = X.text, axis.title = a.text, title = t.text ) +
  theme(strip.background = element_rect(colour = "grey", fill = "white", size = 2, linetype = 1)) +
  labs(x = NULL, y = NULL, title = NULL)  +  
  facet_wrap(~Drug_name1, nrow = 1)



cc +  scale_y_continuous(limits=c(-0,1.2), breaks=c(0, 0.2, 0.4,
                                                    0.6, .8, 1), expand = c(0,0)) +
  scale_x_continuous(limits=c(-5,12), breaks=c(-4,-2,0, 2, 4,
                                               6, 8, 10), expand = c(0,0)) 

rm(dat.sub)



######
# grp2
######

in.y <- which(dat.sub$Drug_name %in% grp2)
dat.sub1 = dat.sub[in.y,]


dat.sub1$Drug_name1 = grp2
dat.c$Drug_name1 = grp2



##########
# Plot parameters
##########
a.text = element_text(size=25, color = "blue")
t.text = element_text(size=25, color = "blue")
Y.text <- element_text( color = "black", size = 19, angle = 0)
X.text <- element_text( color = "black", size = 19, angle = 0)
strip.text.a = element_text(size = 24, color = "black")

dat.sub1$Metadata_dil_factor1 = as.character(dat.sub1$Metadata_dil_factor)


cc = ggplot(dat.sub1, aes(ObjectTotalInten_NUC_A, colour = Metadata_dil_factor1)) +
     geom_density(aes(linetype = Treatment, size = Treatment)) +
     geom_density(data = dat.c, aes(ObjectTotalInten_NUC_A, size = Treatment), colour = "black", linetype = "dashed") +
  
  ##########################
  # mapped variables - linetype, size and color
  scale_linetype_manual(values=c("solid", "solid")) +
  scale_size_manual( values = c(.8,.9) ) +
  ##########################
  #colors
  scale_color_brewer(palette = "Purples") + 
  ##########################

  theme_minimal() +
  theme(legend.position = 'none') +
  theme(strip.text = strip.text.a, axis.text.y = Y.text, 
        axis.text.x = X.text, axis.title = a.text, title = t.text ) +
  theme(strip.background = element_rect(colour = "grey", fill = "white", size = 2, linetype = 1)) +
  labs(x = NULL, y = NULL, title = NULL)  +  
  facet_wrap(~Drug_name1, nrow = 1)

cc +  scale_y_continuous(limits=c(-0,1.2), breaks=c(0, 0.2, 0.4,
                                                    0.6, .8, 1), expand = c(0,0)) +
  scale_x_continuous(limits=c(-5,12), breaks=c(-4,-2,0, 2, 4,
                                               6, 8, 10), expand = c(0,0)) 

rm(dat.sub)


######
# grp3
######

in.y <- which(dat.sub$Drug_name %in% grp3)
dat.sub1 = dat.sub[in.y,]


dat.sub1$Drug_name1 = grp3
dat.c$Drug_name1 = grp3



##########
# Plot parameters
##########
a.text = element_text(size=25, color = "blue")
t.text = element_text(size=25, color = "blue")
Y.text <- element_text( color = "black", size = 19, angle = 0)
X.text <- element_text( color = "black", size = 19, angle = 0)
strip.text.a = element_text(size = 24, color = "black")

dat.sub1$Metadata_dil_factor1 = as.character(dat.sub1$Metadata_dil_factor)


cc = ggplot(dat.sub1, aes(ObjectTotalInten_NUC_A, colour = Metadata_dil_factor1)) +
     geom_density(aes(linetype = Treatment, size = Treatment)) +
     geom_density(data = dat.c, aes(ObjectTotalInten_NUC_A, size = Treatment), colour = "black", linetype = "dashed") +
  
  ##########################
  # mapped variables - linetype, size and color
  scale_linetype_manual(values=c("solid", "solid")) +
  scale_size_manual( values = c(.8,.9) ) +
  ##########################
  #colors
  scale_color_brewer(palette = "Blues") + 
  ##########################
  theme_minimal() +
  theme(legend.position = 'none') +
  theme(strip.text = strip.text.a, axis.text.y = Y.text, 
        axis.text.x = X.text, axis.title = a.text, title = t.text ) +
  theme(strip.background = element_rect(colour = "grey", fill = "white", size = 2, linetype = 1)) +
  
  labs(x = NULL, y = NULL, title = NULL)  +  
  facet_wrap(~Drug_name1, nrow = 1)


cc +  scale_y_continuous(limits=c(-0,1.2), breaks=c(0, 0.2, 0.4,
                                                    0.6, .8, 1), expand = c(0,0)) +
  scale_x_continuous(limits=c(-5,12), breaks=c(-4,-2,0, 2, 4,
                                               6, 8, 10), expand = c(0,0)) 

rm(dat.sub)

######
# grp4
######

in.y <- which(dat.sub$Drug_name %in% grp4)
dat.sub1 = dat.sub[in.y,]


dat.sub1$Drug_name1 = grp4
dat.c$Drug_name1 = grp4

# Due to strong cell loss (Percent of control = 4%) this data is removed for the visualization
ind.t = which(dat.sub1$drug_conc == "tanespimycin_5_1")
# dat.sub1[ind.t,12] = NA
dat.sub1 = dat.sub1[-ind.t,]


##########
# Plot parameters
##########
a.text = element_text(size=25, color = "blue")
t.text = element_text(size=25, color = "blue")
Y.text <- element_text( color = "black", size = 19, angle = 0)
X.text <- element_text( color = "black", size = 19, angle = 0)
strip.text.a = element_text(size = 24, color = "black")

dat.sub1$Metadata_dil_factor1 = as.character(dat.sub1$Metadata_dil_factor)


cc = ggplot(dat.sub1, aes(ObjectTotalInten_NUC_A, colour = Metadata_dil_factor1)) +
     geom_density(aes(linetype = Treatment, size = Treatment)) +
     geom_density(data = dat.c, aes(ObjectTotalInten_NUC_A, size = Treatment), colour = "black", linetype = "dashed") +
  
    ##########################
    # mapped variables - linetype, size and color
    scale_linetype_manual(values=c("solid", "solid")) +
    scale_size_manual( values = c(.8,.9) ) +
    ##########################
    # colors
    # Reds
    scale_colour_manual(values = c("#FEE5D9", "#FCBBA1", "#FC9272", "#FB6A4A", "#EF3B2C", "#CB181D")) +
  
    ##########################
    theme_minimal() +
    theme(legend.position = 'none') +
    theme(strip.text = strip.text.a, axis.text.y = Y.text, 
        axis.text.x = X.text, axis.title = a.text, title = t.text ) +
    theme(strip.background = element_rect(colour = "grey", fill = "white", size = 2, linetype = 1)) +
  
    labs(x = NULL, y = NULL, title = NULL)  +  
    facet_wrap(~Drug_name1, nrow = 1)



cc + scale_y_continuous(limits=c(-0,1.2), breaks=c(0, 0.2, 0.4,
                                                    0.6, .8, 1), expand = c(0,0)) +
  scale_x_continuous(limits=c(-5,12), breaks=c(-4,-2,0, 2, 4,
                                               6, 8, 10), expand = c(0,0)) 

rm(dat.sub)

