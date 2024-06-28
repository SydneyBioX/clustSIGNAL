# clustSIGNAL

An R package to perform spatial clustering on spatially-resolved transcritomics datasets.

clustSIGNAL: clustering of Spatially Informed Gene expression with Neighbourhood Adapted Learning. Here, we calculate entropy as a measure of "domainness" of cell neighbourhoods, and use it to generate weight distributions to perform adaptive smoothing of gene expression. Homogeneous neighbourhoods have low entropy, and so, smoothing is performed over more cells in these neighbourhoods, whereas heterogeneous neighbourhoods have high entropy and are smoothed less. This approach not only overcomes data sparsity in the gene expression data but also incorporates spatial context in the form of cell arrangement information from the neighbourhoods. The resulting adaptively smoothed gene expression data is used for clustering and could be used for other downstream analyses.

## Installation

You can install the clustSIGNAL from
[GitHub](https://github.com/SydneyBioX/clustSIGNAL) with:

``` r
# install.packages("devtools")
devtools::install_github("SydneyBioX/clustSIGNAL")
```
