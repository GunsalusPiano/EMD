

  ########################
  # R script for figure 2g: Assay and experimental design
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
  cell.dat <- read_csv("data_fig2G.csv")
  
  
  # extract 
  unique(cell.dat$drug_conc)
  treats = c("irinotecan_5_0.25", "monensin_0.3125_0.015625", "sirolimus_10_0.5", "vincristine_10_0.5", "DMSO_NULL_NULL")
  t.keep = which(cell.dat$drug_conc %in% treats)
  df.sub1 = cell.dat[t.keep,] #just treatments


  ########################
  # add annotations

  unique(df.sub1$drug_conc)
  unique(df.sub1$Drug_name)
  df.sub1$treatment_label = df.sub1$drug_conc

  indA = which("irinotecan_5_0.25" == df.sub1$drug_conc)
  df.sub1$treatment_label[indA] = "Treatment A"
  indB = which("monensin_0.3125_0.015625" == df.sub1$drug_conc)
  df.sub1$treatment_label[indB] = "Treatment B"
  indC = which("sirolimus_10_0.5" == df.sub1$drug_conc)
  df.sub1$treatment_label[indC] = "Treatment C"
  indD = which("vincristine_10_0.5" == df.sub1$drug_conc)
  df.sub1$treatment_label[indD] = "Treatment D"
  indE = which("DMSO_NULL_NULL" == df.sub1$drug_conc)
  df.sub1$treatment_label[indE] = "Control"




########################
# ggplot parameters

a.text = element_text(size=16, color = "black")
t.text = element_text(size=16, color = "black")
Y.text <- element_text( color = "black", size = 13, angle = 0)
X.text <- element_text( color = "black", size = 13, angle = 0)
strip.text.a = element_text(size = 14, color = "black")



########################
# Plot cell cycle feature
# Density ridge using quartile option

ggplot(df.sub1, aes(x = ObjectTotalInten_NUC_A, y = treatment_label, fill = factor(stat(quantile)))) +
  stat_density_ridges(
    geom = "density_ridges_gradient", calc_ecdf = TRUE,
    quantiles = 4, quantile_lines = TRUE
  ) +
  theme_classic2() +
  scale_fill_viridis_d(name = "Quartiles") +         
  labs(title = NULL, y = NULL, x = "Total Nucleus Intensity (Hoechst) ")





#########################
# Adding annotations for another figure with modified labels

df.sub1$treatment_label1 = df.sub1$drug_conc

indA = which("irinotecan_5_0.25" == df.sub1$drug_conc)
df.sub1$treatment_label1[indA] = "A"

indB = which("monensin_0.3125_0.015625" == df.sub1$drug_conc)
df.sub1$treatment_label1[indB] = "B"

indC = which("sirolimus_10_0.5" == df.sub1$drug_conc)
df.sub1$treatment_label1[indC] = "C"

indD = which("vincristine_10_0.5" == df.sub1$drug_conc)
df.sub1$treatment_label1[indD] = "D"

indE = which("DMSO_NULL_NULL" == df.sub1$drug_conc)
df.sub1$treatment_label1[indE] = "Control"





unique(df.sub1$treatment_label1)
# "Control" "A"       "C"       "D"       "B"      

df=df.sub1
df$method = df$treatment_label1

# ordering the samples 
df$method <- factor(df$method, levels = c("D", "C", "B", "A","Control"), ordered = TRUE)




#########################
# Plot parameters

Y.text <- element_text( color = "black", size = 12, angle = 0)
X.text <- element_text( color = "black", size = 12, angle = 0)
a.text = element_text(size = 15, color = "black")
t.text = element_text(size = 15, color = "black")


########################
# Plot cell cycle feature
# Density ridge using quartile option

ggplot(df, aes(x=ObjectTotalInten_NUC_A, y=method, fill = factor(stat(quantile)))) +
  stat_density_ridges(
    geom = "density_ridges_gradient", calc_ecdf = TRUE,
    quantiles = 4, quantile_lines = TRUE
  ) +
  theme_classic2() +
  theme(legend.position ='top') + 
  
  scale_fill_viridis_d(name = "Quartiles") +         #+ facet_wrap(~Metadata_dil_factor, nrow=1)
  # labs(title = 'Cell cycle feature distributions', y = NULL, x = "Total Nucleus Intensity Hoechst Channel")
  labs(title = NULL, y = NULL, x = "Total Nucleus Intensity") +
  theme(strip.text = strip.text.a, axis.text.y = Y.text, axis.text.x = X.text, axis.title = a.text,
        title = t.text, plot.title = element_text(size=22) ) +
  scale_x_continuous(breaks=seq(-6, 12, 2))


