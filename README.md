## Flow cytometry data analysis
The goal for this repository is to compile a basic data analysis workflow of flow cytometry data that can be used by others in our research lab, or by anyone else, for:
1. anomalous event detection and removal
2. batch effect correction
3. clustering

The main libraries, packages that are used in the current workflow are as follows:
- [Tidyverse](https://www.tidyverse.org)
- [FlowCore](https://bioconductor.org/packages/release/bioc/html/flowCore.html) for handling of flow cytometry data,
  [link for paper](https://doi.org/10.1186/1471-2105-10-106)
- [FlowAI](https://www.bioconductor.org/packages/release/bioc/html/flowAI.html) for quality control, [link for paper](https://doi.org/10.1093/bioinformatics/btw191)
- [cyCombine](https://github.com/biosurf/cyCombine) for batch effect correction, [link for paper](https://doi.org/10.1038/s41467-022-29383-5)
- [FlowSOM](https://www.bioconductor.org/packages/release/bioc/html/FlowSOM.html) for clustering, [link for paper](https://doi.org/10.1002/cyto.a.22625)

The documentation, vignettes and other information of these packages can be found on the above links. Other packages, such as [umap](https://cran.r-project.org/web/packages/umap/index.html), [patchwork](https://cran.r-project.org/web/packages/patchwork/index.html) and [Biobase](https://bioconductor.org/packages/release/bioc/html/Biobase.html), are also used by some specific functions.
<br>
<br>
### Useful information

For learning flow cytometry data handling, cleaning, compensation, gating, transformation and basic plotting in R, I recommend this [tutorial](https://github.com/hally166/R_flowcytometry_course) from Christopher Hall and the following [page](https://med.virginia.edu/flow-cytometry-facility/resources/r-script/).

Marker expression data needs to be transformed before batch effect correction for which the asinh, or inverse hyperbolic sine, transformation is applied in this workflow. Setting an appropriate value of cofactor for the transformation is a crucial step. The following sources may help to figure out what cofactor to use: 
[1.](https://cytoforum.stanford.edu/viewtopic.php?f=3&t=1498) 
[2.](https://github.com/maxentile/advanced-ml-project/issues/2) 
[3.](https://www.researchgate.net/figure/Selecting-the-optimal-value-of-cofactor-using-flowScape-The-distributions-of-CompControl_fig2_224915947)

cyCombine requires the data as a tibble to work with, therefore .fcs files or flowSets have to be converted for which some additional metadata files (metadata_cC.xlsx and panel_cC.xlsx) are needed. Examples of how to prepare these files are shown in the *metadata_preparation_for_cyCombine.pdf*. For more information please see the [Reference manual](https://biosurf.org/cyCombine_ref_manual.html) of cyCombine.

To convert the output of batch effect correction into a flowFrame with the *save_as_ff* function, or into a .fcs file with the *save_as_fcs* function, I utilized some code from Yann Abraham: [link to the original](https://gist.github.com/yannabraham/c1f9de9b23fb94105ca5). Alternatively, the output could be saved as a .csv file and then converted into a .fcs file for which [this](https://floreada.io/fcscreate) might be a useful tool.
