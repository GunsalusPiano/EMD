########################
# function for connecting umap compounds by their concentration
# DIM 1 and DIM 2
########################

make_line_segments <- function(Features_umap_df){
  treatment <- character(0)
  dil_factor <- numeric(0)
  x_start <- numeric(0)
  x_end   <- numeric(0)
  y_start <- numeric(0)
  y_end   <- numeric(0)
  
  for( ii in levels(factor(Features_umap_df$Drug_name)) ){
    temp <- subset(Features_umap_df,Drug_name == ii)
    if(dim(temp)[1] < 2){ next }
    temp <- temp[order(temp$Metadata_dil_factor),]
    
    treatment  <- append(treatment,rep(ii, length( temp$V1[ 1:length(temp$V1)-1 ] ) ) ) 
    dil_factor <- append(dil_factor, temp$Metadata_dil_factor[1:length(temp$V1)-1] )
    x_start <- append(x_start , temp$V1[1:length(temp$V1)-1] )
    x_end   <- append(x_end   , temp$V1[2:length(temp$V1)]   )
    y_start <- append(y_start , temp$V2[1:length(temp$V2)-1] )
    y_end   <- append(y_end   , temp$V2[2:length(temp$V2)]   )
    
  }
  
  umap.lines = data.frame(Drug_name = treatment, Metadata_dil_factor = dil_factor, x1 = x_start, x2 = x_end, y1 = y_start, y2 = y_end)
  
  return(umap.lines)
}






########################
# function for connecting umap compounds by their concentration
# DIM 2 and DIM 3
########################

make_line_segments1 <- function(Features_umap_df){
  treatment <- character(0)
  dil_factor <- numeric(0)
  x_start <- numeric(0)
  x_end   <- numeric(0)
  y_start <- numeric(0)
  y_end   <- numeric(0)
  
  for( ii in levels(factor(Features_umap_df$Drug_name)) ){
    temp <- subset(Features_umap_df,Drug_name == ii)
    if(dim(temp)[1] < 2){ next }
    temp <- temp[order(temp$Metadata_dil_factor),]
    
    treatment  <- append(treatment,rep(ii, length( temp$V2[ 1:length(temp$V2)-1 ] ) ) ) 
    dil_factor <- append(dil_factor, temp$Metadata_dil_factor[1:length(temp$V2)-1] )
    x_start <- append(x_start , temp$V2[1:length(temp$V2)-1] )
    x_end   <- append(x_end   , temp$V2[2:length(temp$V2)]   )
    y_start <- append(y_start , temp$V3[1:length(temp$V3)-1] )
    y_end   <- append(y_end   , temp$V3[2:length(temp$V3)]   )
    
  }
  
  umap.lines1 = data.frame(Drug_name = treatment, Metadata_dil_factor = dil_factor, x1 = x_start, x2 = x_end, y1 = y_start, y2 = y_end)
  
  return(umap.lines1)
}


