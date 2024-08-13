###################################################################
# Function to create UMAPs for batch effect correction evaluation #
###################################################################

# this function uses the uncorrected and corrected flowsets
# if n_samples is drastically decreased might consider increasing marker_alpha for better visibility of the plots

plot_batches_UMAP <- function(uncorr, corr, markers, 
                              n_samples = 50000, 
                              marker_alpha = 0.15, 
                              output_file = "Batch_correction_UMAPs", 
                              grid = F) {
  library(tidyverse)
  library(umap)
  library(patchwork)
  
  # down sample the data
  downSamp <- sample(1:nrow(uncorr), n_samples)
  uncorr_ds <- uncorr[downSamp, ]
  corr_ds <- corr[downSamp, ]
  print("Down sampling is done.")
  print("UMAPs are being made. It might take some time...")
  
  # UMAPs
  uncorr_umap <- umap(uncorr_ds[, markers])
  print("UMAP for the uncorrected data is done...")
  corr_umap <- umap(corr_ds[, markers])
  print("All UMAPs are done.")
  print("Plots are being created...")
  
  # create dataframes for plotting
  uncorr_df <- data.frame(uncorr_umap$layout)
  uncorr_df <- cbind(uncorr_df, batch = uncorr_ds$batch)
  
  corr_df <- data.frame(corr_umap$layout)
  corr_df <- cbind(corr_df, batch = corr_ds$batch)
  
  #plots
  if (grid == T) {
    gridAttr = element_line(color = "lightgrey", linewidth = 0.1)
    axTextAttr = element_text(size = 4)
    axTicksAttr = element_line(color = "lightgrey", linewidth = 0.1)
  } else {
    gridAttr = element_blank()
    axTextAttr = element_blank()
    axTicksAttr = element_blank()
  }
  
  p1 <- uncorr_df %>% 
    ggplot(aes(x = uncorr_df[,1], y = uncorr_df[,2], color = uncorr_df[,3])) +
    geom_point(alpha = marker_alpha, shape = 1, size = 0.05) +
    scale_color_brewer(palette = "Set1") +
    labs(title = "Uncorrected",
         x = "UMAP1",
         y = "UMAP2",
         color = "Batch") +
    theme_light() +
    theme(plot.title = element_text(hjust = 0.5, size = 10),
          axis.title = element_text(size = 6),
          axis.text = axTextAttr,
          axis.ticks = axTicksAttr,
          legend.title = element_text(size = 6),
          legend.text = element_text(size = 4),
          panel.border = element_rect(linewidth = 0.5, color = "black"),
          panel.grid.major = gridAttr,
          panel.grid.minor = element_blank(),
          aspect.ratio = 1,
          plot.margin = unit(c(0,30,0,0), "pt"))
  
  p2 <- corr_df %>% 
    ggplot(aes(x =corr_df[,1], y = corr_df[,2], color = corr_df[,3])) +
    geom_point(alpha = marker_alpha, shape = 1, size = 0.05) +
    scale_color_brewer(palette = "Set1") +
    labs(title = "Corrected",
         x = "UMAP1",
         y = "UMAP2",
         color = "Batch") +
    theme_light() +
    theme(plot.title = element_text(hjust = 0.5, size = 10),
          axis.title = element_text(size = 6),
          axis.text = axTextAttr,
          axis.ticks = axTicksAttr,
          legend.title = element_text(size = 6),
          legend.text = element_text(size = 4),
          panel.border = element_rect(linewidth = 0.5, color = "black"),
          panel.grid.major = gridAttr,
          panel.grid.minor = element_blank(),
          aspect.ratio = 1,
          plot.margin = unit(c(0,0,0,30), "pt"))
  
  p1 | p2
  
  # save the plots
  ggsave(paste0(output_file, ".jpg"), dpi = 600, width = 7, height = 3)
  print("Plots are saved.")
  print("DONE!")
}


############################################################################
# Function to convert the result of batch effect correction to a flowFrame #
############################################################################

# to save batches as separate flowFrames set separate_batches to True

save_as_ff <- function(corr, markers, separate_batches = F) {
  library(Biobase)
  library(flowCore)
  
  if (separate_batches == T) {
    library(tidyverse)
    print("flowFrames are creating...")
    
    batchIDs <- unique(corr$batch) # to get batch IDs
    for (n in batchIDs) {
      # select the corrected data for the actual batch and the markers of interest
      dta <- corr %>% 
        dplyr::filter(batch == n)
      dta <- dta[, markers]
      
      # prepare metadata (required for creating a flowFrame)
      meta <- data.frame(name=dimnames(dta)[[2]],
                         desc=paste(dimnames(dta)[[2]],'marker'))
      meta$range <- apply(apply(dta,2,range),2,diff)
      meta$minRange <- apply(dta,2,min)
      meta$maxRange <- apply(dta,2,max)
      
      # create a flowFrame with the data and metadata
      ff_temp <- new("flowFrame",
                     exprs=data.matrix(dta),
                     parameters=AnnotatedDataFrame(meta))
      assign(paste0("ff_corrected_batch", n), ff_temp, envir = .GlobalEnv) 
    }
  } else {
    print("flowFrame is creating...")
    # select the corrected data for the markers of interest
    dta <- corr[, markers]
    
    # prepare metadata
    meta <- data.frame(name=dimnames(dta)[[2]],
                       desc=paste(dimnames(dta)[[2]],'marker'))
    meta$range <- apply(apply(dta,2,range),2,diff)
    meta$minRange <- apply(dta,2,min)
    meta$maxRange <- apply(dta,2,max)
    
    # create a flowFrame with the data and metadata
    ff_corrected <<- new("flowFrame",
                         exprs=data.matrix(dta),
                         parameters=AnnotatedDataFrame(meta))
  }
  print("DONE!")
}


############################################################################
# Function to convert the result of batch effect correction to a .fcs file #
############################################################################

# to save batches as separate .fcs files set separate_batches to True

save_as_fcs <- function(corr, markers, separate_batches = F) {
  library(Biobase)
  library(flowCore)
  
  if (separate_batches == T) {
    library(tidyverse)
    print("FCS files are creating...")
    
    batchIDs <- unique(corr$batch) # to get batch IDs
    for (n in batchIDs) {
      # select the corrected data for the actual batch and the markers of interest
      dta <- corr %>% 
        dplyr::filter(batch == n)
      dta <- dta[, markers]
      
      # prepare metadata (required for creating a flowFrame)
      meta <- data.frame(name=dimnames(dta)[[2]],
                         desc=paste(dimnames(dta)[[2]],'marker'))
      meta$range <- apply(apply(dta,2,range),2,diff)
      meta$minRange <- apply(dta,2,min)
      meta$maxRange <- apply(dta,2,max)
      
      # create a flowFrame with the data and metadata
      ff_temp <- new("flowFrame",
                     exprs=data.matrix(dta),
                     parameters=AnnotatedDataFrame(meta))
      
      # save flowFrame as a .fcs file
      write.FCS(ff_temp, paste0("corrected_batch", n, ".fcs"))
    }
  } else {
    print("FCS file is creating...")
    # select the corrected data for the markers of interest
    dta <- corr[, markers]
    
    # prepare metadata
    meta <- data.frame(name=dimnames(dta)[[2]],
                       desc=paste(dimnames(dta)[[2]],'marker'))
    meta$range <- apply(apply(dta,2,range),2,diff)
    meta$minRange <- apply(dta,2,min)
    meta$maxRange <- apply(dta,2,max)
    
    # create a flowFrame with the data and metadata
    ff_temp <- new("flowFrame",
                   exprs=data.matrix(dta),
                   parameters=AnnotatedDataFrame(meta))
    
    # save flowFrame as a .fcs file
    write.FCS(ff_temp, "corrected.fcs")
  }
  print("DONE!")
}


#########################################################
# Function to create essential plots of FlowSOM results #
#########################################################

# the res parameter sets the resolution for the saved files in dpi

plot_fSOM_results <- function(fSOM, markers, res = 300) {
  library(tidyverse)
  library(FlowSOM)
  
  # plot the Minimal Spanning Tree (MST)
  PlotStars(fSOM,
            backgroundValues = fSOM$metaclustering,
            equalNodeSize = T,
            view = "MST")
  ggsave("flowSOMresults_MST.jpg", dpi = res)
  
  # plot the MST with representative node sizes
  PlotStars(fSOM,
            backgroundValues = fSOM$metaclustering,
            equalNodeSize = F,
            view = "MST")
  ggsave("flowSOMresults_MST2.jpg", dpi = res)
  
  # plot the grid representation
  PlotStars(fSOM,
            backgroundValues = fSOM$metaclustering,
            equalNodeSize = F,
            view = "grid")
  ggsave("flowSOMresults_grid.jpg", dpi = res)
  
  # labels for the MST
  PlotLabels(fSOM, labels = fSOM$metaclustering)
  ggsave("flowSOMresults_clusterLabels.jpg", dpi = res)
  
  # plot the markers on the MST
  for (m in markers) {
    PlotMarker(fSOM, m)
    ggsave(paste0("flowSOMresults_MST_", m, ".jpg"), dpi = res)
  }
  
  # get and reorganize data for heatmap
  htm_data <- GetMetaclusterMFIs(fSOM)
  htm_data <- round(htm_data, digits = 2)
  htm_data <- rownames_to_column(htm_data, var = "cluster")
  
  # plot heatmap
  htm_data %>% 
    pivot_longer(cols = !cluster, names_to = "marker", values_to = "MFI") %>% 
    ggplot(aes(x = marker, y = cluster, fill = MFI, label = MFI)) +
    geom_tile() +
    geom_text() +
    scale_fill_gradientn(colours = c("darkblue", "blue", "cyan1", "springgreen",
                                     "yellow", "red", "darkred")) +
    labs(x = "Marker",
         y = "Cluster",
         title = "Median fluorescence intensity") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid = element_blank(),
          aspect.ratio = 1)
  ggsave("flowSOMresults_heatmap.jpg", dpi = res)
  
  # get cell number per cluster
  nCell <- as.data.frame(GetCounts(fSOM))
  nCell <- rownames_to_column(nCell, var = "cluster")
  
  # plot cell number per cluster
  nCell %>%
    ggplot(aes(x = cluster, y = GetCounts(fSOM))) +
    geom_col(fill = "#79dbc0") +
    labs(x = "Cluster",
         y = "Count",
         title = "Number of cells per (meta)cluster") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "NULL",
          aspect.ratio = 1)
  ggsave("flowSOMresults_cellCount.jpg", dpi = res)
  
  print("DONE!")
}
