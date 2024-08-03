library(tidyverse)
library(flowCore)
library(flowAI)
library(cyCombine)
library(FlowSOM)

# The data used here was preprocessed in FlowJo: ...


################################
# Work with a single .fcs file #
################################

# import a .fcs file into a flowframe (ff) to work with
# place file in the working directory and define file name in double quotes 
ff <- flowCore::read.FCS("STATIC_FULL_1.fcs",
                         transformation = F)

# several parameters of a ff can be explored
ff
summary(ff)
parameters(ff)
colnames(ff)
names(ff)
featureNames(ff)
head(exprs(ff))
each_col(ff, median)
keyword(ff)

# clean the data with flowAI
# flowAI does QC along 3 parameters: flow rate, signal acquisition and dynamic range
# by default, it uses all 3 parameters but it can be specified otherwise (see manual)
# output files with QC reports are saved in a new directory, named resultsQC
ff_clean <- flow_auto_qc(ff)

# there is also a Shiny application to perform interactive QC with flowAI
# run the code below to start the app
flow_iQC()


################################
# Work with several .fcs files #
################################

# .fcs files can be imported into one flowset (fs) to make handling easier
# define file path (in double quotes) to the directory where files are stored
file_path = "C:/Users/peti/Documents/flow_cytometry_workflow"
files <- list.files(path = file_path, pattern = ".fcs$")
fs <- read.flowSet(files, path = file_path, 
                   transformation = F)

# content of a fs
fs
summary(fs)
sampleNames(fs)
fsApply(fs, colnames)
fsApply(fs, exprs)
fsApply(fs, each_col, median)
fsApply(fs, keyword)
# a single element of a fs can be selected and can be explored as below
# select an element of a fs by indexing
fs[[1]]
parameters(fs[[1]])
colnames(fs[[1]])
names(fs[[1]])
featureNames(fs[[1]])
head(exprs(fs[[1]]))
each_col(fs[[1]], median)
keyword(fs[[1]])

# clean the fs with flowAI
fs_clean <- flow_auto_qc(fs)


##########################################
# Batch effect correction with cyCombine #
##########################################

# the above created cleaned fs can be used as an input for cyCombine

# alternatively, other .fcs files can be read into a fs to work with
fs_cC <- compile_fcs(data_dir = file_path,
                     pattern = "\\.fcs")

# cyCombine requires the fs to be converted into a tibble
# in addition, metadata and panel files are needed for this step
# (examples of these files can be found attached)
df_cC <- convert_flowset(flowset = fs_clean, 
                         metadata = file.path(file_path, "metadata_cC.xlsx"), 
                         sample_ids = NULL, 
                         batch_ids = "batch", 
                         filename_col = "filename", 
                         condition = "condition", 
                         down_sample = F,
                         panel = file.path(file_path, "panel_cC.xlsx"), 
                         panel_channel = "channel",
                         panel_antigen = "antigen")

# define markers of interest
# the rest will be ignored for transformation and for batch effect correction
markers <- c("CD27", "IgD", "CD38")

# transform the data for batch effect correction
# for conventional flow cytometry data, cofactor is usually set to 150-250
uncorrected <- transform_asinh(df_cC, 
                               markers = markers,
                               cofactor = 150,
                               derand = F)

# store the results before correction (optional)
saveRDS(uncorrected, file = file.path(file_path, "uncorrected.RDS"))

# check for batch effect prior to correction (optional)
detect_batch_effect(uncorrected,
                    batch_col = 'batch',
                    out_dir = paste0(file_path, "/batch_effect_check"),
                    xdim = 7,
                    ydim = 7,
                    seed = 434)

# create cell type labels by using a self-organizing map (SOM)
# the dimensions of the SOM grid are defined by xdim and ydim
# a larger grid creates more clusters and a finer separation of cell groups
# rlen shows the number of times the data is presented to the SOM network
labels <- uncorrected %>% 
  normalize(markers = markers,
            norm_method = "scale",   # can also be "rank" for stronger correction
            ties.method = "average") %>%   # can also be "minimum"
  create_som(markers = markers,
             rlen = 10,   # higher values may lead to better results
             seed = 101,
             xdim = 7,
             ydim = 7)

# run batch correction
corrected <- uncorrected %>%
  correct_data(label = labels,
               covar = "condition",
               markers = markers)

# store the results after correction (optional)
saveRDS(corrected, file.path(file_path, "corrected.RDS"))

# visualize the results of batch correction
plot_batches_UMAP(uncorrected, corrected, markers)

plot_density(uncorrected,
             corrected,
             markers = markers,
             filename = "batch_correction_densityPlots.jpg",
             y = "batch",
             ncol = 3)   # should be equal to the number of markers

# convert the results of batch correction into a ff
# to save batches as separate ff-s set separate_batches to True
# output: ff_corrected or ff_corrected_batchx/y/z
save_as_ff(corrected, markers)
save_as_ff(corrected, markers, separate_batches = T)

# alternatively, results can be saved as .fcs files as well
# output: corrected.fcs or corrected_batchx/y/z.fcs
save_as_fcs(corrected, markers)
save_as_fcs(corrected, markers, separate_batches = T)


###########################
# Clustering with flowSOM #
###########################

# ff-s of the batch corrected data (e.g. ff_corrected) can be used with flowSOM
# create a flowSOM object (fSOM) which contains the results of clustering
# the dimensions of the SOM grid are defined by xdim and ydim
fSOM <- FlowSOM(ff_corrected,
                compensate = F,
                transform = F,
                scale = F,   # define any if further scaling is necessary
                colsToUse = NULL,   # set to NULL to use all columns for clustering
                rlen = 10,   # higher values may lead to better results
                xdim = 8,
                ydim = 8,
                nClus = 6,   # number of metaclusters
                #maxMeta = NULL,   # if no exact number of metaclusters is required
                seed = 104)

# explore the results of clustering
print(fSOM)

GetClusters(fSOM)   # cluster labels
GetClusterCVs(fSOM)   # coefficient for variation of each cluster
GetClusterMFIs(fSOM)   # median values of each cluster

GetMetaclusters(fSOM)
GetMetaclusterCVs(fSOM)
GetMetaclusterMFIs(fSOM)

# plot a summary of the fSOM object
FlowSOMmary(fSOM)

# a simple way to create essential plots of clustering results
plot_fSOM_results(fSOM, markers)