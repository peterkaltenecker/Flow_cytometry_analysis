# Function to create UMAPs for comparison
# This function uses the uncorrected and corrected flowsets
# If n_samples is drastically decreased might consider increasing marker_alpha for better visibility of the plots

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
  uncorr.umap <- umap(uncorr_ds[, markers])
  print("UMAP for the uncorrected data is done...")
  corr.umap <- umap(corr_ds[, markers])
  print("All UMAPs are done.")
  print("Plots are being created...")
  
  # create dataframes for plotting
  uncorr_df <- data.frame(uncorr.umap$layout)
  uncorr_df <- cbind(uncorr_df, batch = uncorr_ds$batch)
  
  corr_df <- data.frame(corr.umap$layout)
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
