
  ########################
  # R script for figure 3B: Summary of positional effect detection
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
  library(car)
  library(viridis)


  ########################
  # Import all output from the two-way ANOVA analysis, summarize results in figures   
  ########################

  # anova_A = twoway_anova_panelA.csv
  # anova_B = twoway_anova_panelB.csv
  # anova_C1 = twoway_anova_panelC1.csv
  # anova_C2 = twoway_anova_panelC2.csv


  anova_A$marker = anova_A$Channel
  anova_B$marker = anova_B$Channel
  anova_C1$marker = anova_C1$Channel
  anova_C2$marker = anova_C2$Channel

  ########################
  # Add stain name annotation for panel A
  ########################

  ch = unique(anova_A$Channel)
  ind1 = which("NUC_A" == anova_A$Channel)
  anova_A$marker[ind1] = "Hoechst_panelA"
  ind2 = which("RNA" == anova_A$Channel)
  anova_A$marker[ind2] = "SYTO14"
  ind3 = which("PMG" == anova_A$Channel)
  anova_A$marker[ind3] = "WGA-AlexaFluor555"
  ind4 = which("MITO" == anova_A$Channel)
  anova_A$marker[ind4] = "MitoTrackerDR"
  ind5 = which("NUC_Tex" == anova_A$Channel)
  anova_A$marker[ind5] = "Hoechst"
  ind6 = which("CELL_COUNT_A" == anova_A$Channel)
  anova_A$marker[ind6] = "Cell_count_panelA"

  ########################
  # Add stain name annotation for panel B
  ########################

  ch = unique(anova_B$Channel)
  ind1 = which("NUC_B" == anova_B$Channel)
  anova_B$marker[ind1] = "Hoechst_panelB"
  ind2 = which("Lyso" == anova_B$Channel)
  anova_B$marker[ind2] = "LysoTrackerRedDND99"
  ind3 = which("Perox" == anova_B$Channel)
  anova_B$marker[ind3] = "PeroxisomeGFP"
  ind4 = which("Lipids" == anova_B$Channel)
  anova_B$marker[ind4] = "LipidTOXDR" # HCS LipidTOX™ Deep Red
  ind5 = which("CELL_COUNT_B" == anova_B$Channel)
  anova_B$marker[ind5] = "Cell_count_panelB"

  ########################
  # Add stain name annotation for panel C1
  ########################

  ch = unique(anova_C1$Channel)
  ind1 = which("NUC_C1" == anova_C1$Channel)
  anova_C1$marker[ind1] = "DRAQ5"
  ind2 = which("ER" == anova_C1$Channel)
  anova_C1$marker[ind2] = "ER-TrackerBlueWhiteDPX"
  ind3 = which("CELL_COUNT_C1" == anova_C1$Channel)
  anova_C1$marker[ind3] = "Cell_count_panelC1"

  ########################
  # Add stain name annotation for panel C2
  ########################

  ch = unique(anova_C2$Channel)
  ind1 = which("NUC_C2" == anova_C2$Channel)
  anova_C2$marker[ind1] = "Hoechst_panelC2"
  ind2 = which("TUB" == anova_C2$Channel)
  anova_C2$marker[ind2] = "Tubulin-GFP"
  ind3 = which("ACT" == anova_C2$Channel)
  anova_C2$marker[ind3] = "Actin-RFP"
  ind4 = which("CELL_COUNT_C2" == anova_C2$Channel)
  anova_C2$marker[ind4] = "Cell_count_panelC2"


  ########################
  # check data
  ########################

  length(unique(anova_A$feature)) 
  length(unique(anova_B$feature)) 
  length(unique(anova_C1$feature)) 
  length(unique(anova_C2$feature)) 
  
  ########################
  # Bind the data
  ########################
  
  anova_all = rbind(anova_A, anova_B, anova_C1, anova_C2)

  unique(anova_all$Channel)
  unique(anova_all$marker)
  
  ########################
  # Order the data by p-value
  ########################
  
  dat.all = anova_all 
  dat.ord = dat.all[order(dat.all$P_r),]
  dat.ord$x = 1:1068

    
    ########################
    # Prepare data for visualizations
    ########################
    
  
  ind1 = which(grepl("1A1",dat.all$plate)|grepl("1B1",dat.all$plate)|grepl("1_C1_1",dat.all$plate)|grepl("1_C2_1",dat.all$plate))
  dat.all$plate_rep = dat.all$plate
  dat.all$plate_rep[ind1] = "plate1_rep1"
  
  ind2 = which(grepl("1A2",dat.all$plate)|grepl("1B2",dat.all$plate)|grepl("1_C1_2",dat.all$plate)|grepl("1_C2_2",dat.all$plate))
  dat.all$plate_rep[ind2] = "plate1_rep2"
  
  ind3 = which(grepl("1A3",dat.all$plate)|grepl("1B3",dat.all$plate)|grepl("1_C1_3",dat.all$plate)|grepl("1_C2_3",dat.all$plate))
  dat.all$plate_rep[ind3] = "plate1_rep3"
  
  ind4 = which(grepl("2A1",dat.all$plate)|grepl("2B1",dat.all$plate)|grepl("2_C1_1",dat.all$plate)|grepl("2_C2_1",dat.all$plate))
  dat.all$plate_rep[ind4] = "plate2_rep1"
  
  ind5 = which(grepl("2A2",dat.all$plate)|grepl("2B2",dat.all$plate)|grepl("2_C1_2",dat.all$plate)|grepl("2_C2_2",dat.all$plate))
  dat.all$plate_rep[ind5] = "plate2_rep2"
  
  ind6 = which(grepl("2A3",dat.all$plate)|grepl("2B3",dat.all$plate)|grepl("2_C1_3",dat.all$plate)|grepl("2_C2_3",dat.all$plate))
  dat.all$plate_rep[ind6] = "plate2_rep3"
  
  ########################
  # Add annotations to facilitate the plots
  ########################
  
    dat.all$marker_sub = dat.all$marker
    ind1 = which(grepl("Cell_count", dat.all$marker))
    dat.all$marker_sub[ind1] = "Cell_count"
    ind1 = which(grepl("Hoechst", dat.all$marker))
    dat.all$marker_sub[ind1] = "Hoechst"
    unique(dat.all$marker_sub)
    ind.count = which(grepl("COUNT",dat.all$feature))
    count_check = dat.all[ind.count,]
    dat.all.new = dat.all
    dat.all.new$feature[ind.count] = dat.all$Channel[ind.count]
  

      ########################
      # Visualizations for negative-log(p-value) 
      ########################
  
        dat.all.new$log_neg = -(dat.all.new$P_r_log)
        unique(dat.all.new$plate_rep)
        
        # Adding labels for plate # and replicate #
        p1r1 = which(dat.all.new$plate_rep == "plate1_rep1")
        p1r2 = which(dat.all.new$plate_rep == "plate1_rep2")
        p1r3 = which(dat.all.new$plate_rep == "plate1_rep3")
        p2r1 = which(dat.all.new$plate_rep == "plate2_rep1")
        p2r2 = which(dat.all.new$plate_rep == "plate2_rep2")
        p2r3 = which(dat.all.new$plate_rep == "plate2_rep3")
        
        # initialize new column
        dat.all.new$plate_rep1 = dat.all.new$plate_rep
        
        # editing lables
        dat.all.new$plate_rep1[p1r1] = "p1.r1"
        dat.all.new$plate_rep1[p1r2] = "p1.r2"
        dat.all.new$plate_rep1[p1r3] = "p1.r3"
        
        dat.all.new$plate_rep1[p2r1] = "p2.r1"
        dat.all.new$plate_rep1[p2r2] = "p2.r2"
        dat.all.new$plate_rep1[p2r3] = "p2.r3"
        
        ########################
        # change labels from channel/marker names to cell structure
        ########################
        
        unique(dat.all.new$marker_sub)
        
        ind1 = which(dat.all.new$marker_sub == "Cell_count")
        ind2 = which(dat.all.new$marker_sub == "Hoechst")
        ind3 = which(dat.all.new$marker_sub == "SYTO14")
        
        ind4 = which(dat.all.new$marker_sub == "WGA-AlexaFluor555")
        ind5 = which(dat.all.new$marker_sub == "MitoTrackerDR")
        ind6 = which(dat.all.new$marker_sub == "LysoTrackerRedDND99")
        
        ind7 = which(dat.all.new$marker_sub == "PeroxisomeGFP")
        ind8 = which(dat.all.new$marker_sub == "LipidTOXDR")
        ind9 = which(dat.all.new$marker_sub == "DRAQ5")
        
        ind10 = which(dat.all.new$marker_sub == "ER-TrackerBlueWhiteDPX")
        ind11 = which(dat.all.new$marker_sub == "Tubulin-GFP")
        ind12 = which(dat.all.new$marker_sub == "Actin-RFP")
        
        ########################
        # initialize new column
        dat.all.new$cell_struc = dat.all.new$marker_sub
        
        ########################
        # editing lables
        dat.all.new$cell_struc[ind1] = "Cell count"
        dat.all.new$cell_struc[ind2] = "DNA (Hoechst)"
        dat.all.new$cell_struc[ind3] = "RNA"
        
        dat.all.new$cell_struc[ind4] = "Golgi-membranes"
        dat.all.new$cell_struc[ind5] = "Mitochondria"
        dat.all.new$cell_struc[ind6] = "Lysosomes"
        
        dat.all.new$cell_struc[ind7] = "Peroxisomes"
        dat.all.new$cell_struc[ind8] = "Lipids"
        dat.all.new$cell_struc[ind9] = "DNA (DRAQ5)"
        
        dat.all.new$cell_struc[ind10] = "ER"
        dat.all.new$cell_struc[ind11] = "Tubulin"
        dat.all.new$cell_struc[ind12] = "Actin"
      
        
        
        ########################
        # Use high/low range for plotting
        ########################
        
        rr.c = range(dat.ord$P_c_log)
        rr.r = range(dat.ord$P_r_log)
        rr.rc = range(dat.ord$P_c_log,dat.ord$P_r_log)
        
        low = rr.rc[1] - 3 
        high = rr.rc[2] + 3
        
        
        ########################
        # ggplot parameters
        ########################
        Y.text <- element_text( color = "black", size = 14, angle = 0)
        X.text <- element_text( color = "black", size = 14, angle = 0)
        a.text = element_text(size = 16, color = "black")
        t.text = element_text(size = 16 ,color = "black")
        strip.text.a = element_text(size = 18, color = "black")
        
        
        ########################
        # Figure 3b  Summary of two-way ANOVA test for row effects
        ########################
        
        ggplot(dat.all.new) + 
        aes(x = plate_rep1, y = log_neg, group = feature, colour = feature_type) +
        geom_line(linetype = 1, size = .6) + 
        geom_hline(yintercept=-log(0.00001), linetype = "dashed", color = "black", size = 0.6) +
        geom_hline(yintercept=-log(0.001), linetype = "dotdash", color = "black", size = 0.6) +
        
        ylim(-high, -low) +
        theme_bw() +
        theme(legend.position = "none") +
        labs(x = NULL, y = "-Log(p-value) Row Effects") + 
        theme(strip.background = element_rect(fill="white"), strip.text = strip.text.a, axis.text.x = X.text, axis.text.y = Y.text, axis.title = a.text, title = t.text) +
        scale_color_manual(values = c("#C3D7A4", "#52854C")) + 
        theme(legend.text=element_text(size=18)) +
        facet_wrap(~ cell_struc,nrow=3) 
      
      
      
      
        ########################
        # Plot legend
        ########################
        
        ########################
        # Create empty plot
        plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
        
        ########################
        # Choose colors
        ch.col = c("#C3D7A4", "#52854C")
        
        legend("center", 
               legend = c("Intensity features", "Morphological features"), 
               # col = c(rgb(0.2,0.4,0.1,0.7), 
               #         rgb(0.8,0.4,0.1,0.7)), 
               col = ch.col,
               #pch = c(19,19, 19, 19, 19, 19, 19, 19, 19,19 ), 
               pch = c(15, 15), 
               
               bty = "n", 
               pt.cex = 2, 
               cex = 1.2, 
               text.col = "black", 
               horiz = F , 
               inset = c(0.1, 0.1))
      
      
      
      
      
      
      
      
      
      