

  ########################
  # R script for figure  5 B & C: Replicate aggregation and EMD profiling using global controls.
  ########################


  ########################
  # Clear variables
  # Set working directory
  ########################

  rm(list=ls())
  
  

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


  ########################
  # Import single cell data - feature: Area of the nucleus
  ########################
  

  ########################
  # set working directory
  # setwd("~/Documents/...")
  
  
  new_trt = read_csv("treated_cells_data_fig5B-C.csv")
  new_ctrl = read_csv("control_cells_data_fig5B-C.csv")
  
  ######################## 
  # ggplot figure parameters
  ########################

  Y.text <- element_text( color = "black", size = 20, angle = 0)
  X.text <- element_text( color = "black",  size = 20, angle = 0) 
  a.text = element_text(size=18, color = "black")
  t.text = element_text(size=10, color = "black")
  strip.text.a = element_text(size = 20, color = "black")

  ########################
  # tc/tc1: treatment vs. control CDFs
  # tcd/tcd1 : treatment vs. control density curves
  # cc/cc1: control vs. control CDFs
  # cd/cd1: control vs. control desnity curves
  ########################
  
  ########################
  # treatment vs. control - density curves
  ########################

tcd = ggplot(new_trt, aes(x = ObjectArea_NUC_A , group = plate_well), col = 'lightgrey') +
    geom_density(size = .2, alpha = .4, col = "grey") + 
  
    # Plot parameters
    # ylim(0, 1.4) +
    # ylim(0, 0.86) +
    
    theme_light() + 
    theme(legend.position = 'none', axis.text.y = Y.text, axis.text.x = X.text, axis.title = a.text, 
          title = t.text, plot.title = element_text(size = 22) ) +
    labs(x = "Area of nucleus (Hoechst channel)", title = "Treated cell populations", y = "Density") 

    # Add the contro-DMSO density curve
    tcd1 = tcd + 
    geom_density(data = new_ctrl, aes(x = ObjectArea_NUC_A , group = Treatment), size = 2, col = "black") 
    tcd1

    ########################
    # treatment vs. control - CDF curves
    ########################


tc = ggplot(new_trt, aes(x=ObjectArea_NUC_A , group = drug_conc), col = 'lightgrey') +
    stat_ecdf(geom = "step", size = .2, alpha = .4, col = "grey") + 
  
    # Plot parameters
    theme_light() + 
    theme(legend.position='none', axis.text.y = Y.text, axis.text.x = X.text, 
          axis.title = a.text, title = t.text, plot.title = element_text(size = 22) ) +
    labs(x = "Area of nucleus (Hoechst channel)", title = "Treated cell populations", y = "Empirical CDFs") 
    
    # Add the contro-DMSO CDF curve
    tc1 = tc + 
    stat_ecdf(data = new_ctrl, geom = "step", aes(x=ObjectArea_NUC_A , group = Treatment), size = 1, col = "black") 
    tc1


    ########################
    # Control vs. control - CDF curves
    ########################


cc = ggplot(new_ctrl, aes(x = ObjectArea_NUC_A , group = plate_well), col = 'lightgrey') +
    stat_ecdf(geom = "step", size = .2, alpha = .4, col = "grey") + 
    
    # Plot parameters
    theme_light() + 
    theme(legend.position='none', axis.text.y = Y.text, axis.text.x = X.text, 
          axis.title = a.text, title = t.text, plot.title = element_text(size=22) ) +
    labs(x = "Area of nucleus (Hoechst channel)", title = "Untreated cell populations", y = "Empirical CDFs") 
     
    # Add the contro-DMSO CDF curve
    cc1 = cc + 
    stat_ecdf(data = new_ctrl, geom = "step", aes(x=ObjectArea_NUC_A , group = Treatment), size = 1, col = "black") 
    cc1
    
    
    ########################
    # Control vs. control - Density curves
    ########################
    
cd = ggplot(new_ctrl, aes(x=ObjectArea_NUC_A , group = plate_well), col = 'lightgrey') +
    geom_density(size=.2, alpha=.4, col = "grey") + theme_light() + 

    ylim(0, .86) +
    theme(legend.position='none', axis.text.y = Y.text, axis.text.x = X.text, axis.title = a.text, title = t.text, plot.title = element_text(size=22) ) +
    labs(x = "Area of nucleus (Hoechst channel)", title = "Untreated cell populations", y = "Density") 
    
    # Add the control curve
    cd1 = cd + 
    geom_density(data=new_ctrl, aes(x=ObjectArea_NUC_A , group = Treatment), size = 2, col = "black") 
    cd1

    ########################
    # Plot multiple figures together
    # grid.arrange(t1,a1,ncol = 1)





    ########################
    # Summary of sample sizes
    ########################
    
    library("dplyr")
    
    # DMSO-control counts
    ctrl.count = new_ctrl %>% count(plate_well)
    summary(ctrl.count)
    
    # Counts in each dmso-control well
    counts.c = as.data.frame(ctrl.count)
    # Sum of counts in the control wells used for global control distributions
    sum(counts.c$n) # 265,638
    
    # chemically treated wells + merged replicates
    trt.count = new_trt %>% count(drug_conc)
    summary(trt.count)
  
    # Per treatment counts, here counts reflect three wells combined into one population
    counts.t = as.data.frame(trt.count)
    
    
    
    # setwd("~/Desktop/Glenn_nov16_2019/manuscript_drafts/manuscript_code/data_files_for_code")
    # setwd("~/Desktop/Glenn_nov16_2019/manuscript_drafts/manuscript_code/data_files_for_code")
    
    # Save to per figure folder "Figure5"
    setwd("~/Documents/NYU_work/Sept_17_2018/R_code_b/Manuscript_Comm_Bio_scripts/Figure5")
    write.csv(counts.t, file = "Rawcellcounts_trt_fig5c.csv")
    write.csv(counts.c, file = "Rawcellcounts_dmso_fig5b.csv")
    
    # Save to Data_files folder
    setwd("~/Documents/NYU_work/Sept_17_2018/R_code_b/Manuscript_Comm_Bio_scripts/Data_files/Fig5")
    write.csv(counts.t, file = "Rawcellcounts_trt_fig5c.csv")
    write.csv(counts.c, file = "Rawcellcounts_dmso_fig5b.csv")
    
    
    
    