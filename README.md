## Flow cytometry data analysis
The goal of this repository is to compile a basic data analysis workflow of flow cytometry data that can be used by others in our research lab, or by any other individal, for:
1. anomalous event detection and removal
2. batch effect correction
3. clustering

The main libraries, packages that are used in the current workflow are as follows:
- tidyverse
- FlowCore for handling of flow cytometry data
- FlowAI
- cyCombine for batch effect correction
- [FlowSOM](https://www.bioconductor.org/packages/release/bioc/html/FlowSOM.html) for clustering, [link for paper](https://doi.org/10.1002/cyto.a.22625)
