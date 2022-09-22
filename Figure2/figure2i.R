

  ########################
  # R script for figure 2i: Assay and experimental design
  ########################

  ########################
  # load libraries
  library(easypackages)
  library(ggridges)
  libraries("rmarkdown","resample","assertthat","lattice","gplots","ggplot2","fields","readr","data.table","grid","MASS", "reshape", "plyr", "ggpubr", "tictoc")  
  libraries("plotly","gridExtra","Rmisc","Hmisc","doBy","caret","corrplot","GGally", "easypackages")  




  ########################
  # import single cell data
  cell.dat <- # data_fig2i.csv
  

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




