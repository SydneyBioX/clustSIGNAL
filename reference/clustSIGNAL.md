# ClustSIGNAL

Clustering method for spatially-resolved cell-state classification of
spatial transcriptomics data. The tool generates and uses an adaptively
smoothed, spatially-informed gene expression for clustering.

## Usage

``` r
clustSIGNAL(
  spe,
  samples,
  dimRed_init = "None",
  dimRed_f = c("None", "embed.smooth"),
  batch = FALSE,
  batch_by = "None",
  NN = 30,
  kernel = c("G", "E"),
  spread = 0.3,
  sort = TRUE,
  threads = 1,
  outputs = c("c", "n", "s", "a"),
  clustParams = list(clust_c = 0, subclust_c = 0, iter.max = 30, k = 10, cluster.fun =
    "louvain")
)
```

## Arguments

- spe:

  a SpatialExperiment object containing spatial coordinates in
  'spatialCoords' matrix and normalised gene expression in 'logcounts'
  assay.

- samples:

  a character indicating name of colData(spe) column containing sample
  names.

- dimRed_init:

  a character indicating the name of the reduced dimensions in the
  SpatialExperiment object (i.e., from reducedDimNames(spe)) to use for
  initial clustering step. Default value is 'None'.

- dimRed_f:

  a character indicating the name of the reduced dimensions in the
  SpatialExperiment object (i.e., from reducedDimNames(spe)) to use for
  final clustering step. Two valid options are "None" (default), which
  triggers a PCA run on smoothed expression, and "embed.smooth", which
  triggers a search for externally-generated "embed.smooth" low
  embedding in reducedDimNames(spe).

- batch:

  a logical parameter for whether to perform batch correction. Default
  value is FALSE.

- batch_by:

  a character indicating name of colData(spe) column containing the
  batch names. Default value is 'None'.

- NN:

  an integer for the number of neighbouring cells the function should
  consider. The value must be greater than or equal to 1. Default value
  is 30.

- kernel:

  a character for type of distribution to be used. The two valid values
  are "G" or "E" for Gaussian and exponential distributions,
  respectively. Default value is "G".

- spread:

  a numeric value for distribution spread, represented by standard
  deviation for Gaussian distribution and rate for exponential
  distribution. Default value is 0.3 for Gaussian distribution. The
  recommended value is 5 for exponential distribution.

- sort:

  a logical parameter for whether to sort the neighbourhood by initial
  clusters. Default value is TRUE.

- threads:

  a numeric value for the number of CPU cores to be used for the
  analysis. Default value set to 1.

- outputs:

  a character for the type of output to return to the user. "c" for data
  frame of cell IDs and their respective ClustSIGNAL cluster labels, "n"
  for ClustSIGNAL cluster dataframe plus neighbourhood matrix, "s" for
  ClustSIGNAL cluster dataframe plus final SpatialExperiment object, or
  "a" for all 3 outputs. Default value is 'c'.

- clustParams:

  a list of parameters for TwoStepParam clustering methods: clust_c is
  the number of centers to use for clustering with KmeansParam. By
  default set to 0, in which case the method uses either 3000 centers or
  1/5th of the total cells in the data as the number of centers,
  whichever is lower. subclust_c is the number of centers to use for
  sub-clustering the initial clusters with KmeansParam. The default
  value is 0, in which case the method uses either 1 center or half of
  the total cells in the initial cluster as the number of centers,
  whichever is higher. iter.max is the maximum number of iterations to
  perform during clustering and sub-clustering with KmeansParam. Default
  value is 30. k is a numeric value indicating the k-value used for
  clustering and sub-clustering with NNGraphParam. Default value is 10.
  cluster.fun is a character indicating the graph clustering method used
  with NNGraphParam. By default, the Louvain method is used.

## Value

a list of outputs depending on the type of outputs specified in the main
function call.

1\. clusters: a data frame of cell names and their ClustSIGNAL cluster
classification.

2\. neighbours: a character matrix containing cells IDs of each cell's
NN neighbours.

3\. spe_final: a SpatialExperiment object containing the original spe
object data plus initial cluster and subcluster labels, entropy values,
smoothed gene expression, and ClustSIGNAL cluster labels.

## Examples

``` r
data(ClustSignal_example)

names(colData(spe))
#> [1] "uniqueID"       "sample_id"      "initCluster"    "initSubcluster"
#> [5] "entropy"        "ClustSIGNAL"   
# identify the column name with sample labels
samples = "sample_id"
res_list <- clustSIGNAL(spe, samples, outputs = "c")
#> 06:30:37 ClustSIGNAL running.
#> 06:30:37 Calculating PCA.
#> 06:30:37 Initial clustering performed. Clusters = 3
#> 06:30:38 Initial sub-clustering performed. Subclusters = 7 
#> 06:30:38 Neighbourhoods defined.
#> 06:30:38 Neighbourhood heterogeneity calculated.
#> 06:30:38 Smoothing performed. NN = 30, Kernel = G, Spread = 0.300000
#> 06:30:38 Calculating PCA using smoothed data.
#> 06:30:39 Final clustering performed on smoothed data. Clusters = 4 
#> 06:30:39 ClustSIGNAL completed.
#> Time difference of 1.579595 secs
```
