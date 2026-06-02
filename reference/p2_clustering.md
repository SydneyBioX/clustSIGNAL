# Final non-spatial clustering

A function to perform clustering on adaptively smoothed gene expression
data to generate ClustSIGNAL clusters.

## Usage

``` r
p2_clustering(
  spe,
  dimRed_f = c("None", "embed.smooth"),
  batch = FALSE,
  batch_by = "None",
  threads = 1,
  clustParams = list(clust_c = 0, subclust_c = 0, iter.max = 30, k = 10, cluster.fun =
    "louvain")
)
```

## Arguments

- spe:

  SpatialExperiment object containing the adaptively smoothed gene
  expression.

- dimRed_f:

  a character indicating the name of the reduced dimensions in the
  SpatialExperiment object (i.e., from reducedDimNames(spe)) to use for
  final clustering step. Two valid options are "None" (default), which
  triggers a PCA run on smoothed expression, and "embed.smooth", which
  triggers a search for "embed.smooth" low embedding in
  reducedDimNames(spe).

- batch:

  a logical parameter for whether to perform batch correction. Default
  value is FALSE.

- batch_by:

  a character indicating name of colData(spe) column containing the
  batch names. Default value is 'None'.

- threads:

  a numeric value for the number of CPU cores to be used for the
  analysis. Default value set to 1.

- clustParams:

  a list of parameters for TwoStepParam clustering methods: clust_c is
  the number of centers to use for clustering with KmeansParam. By
  default set to 0, in which case the method uses either 3000 centers or
  1/5th of the total cells in the data as the number of centers,
  whichever is lower. subclust_c is the number of centers to use for
  sub-clustering the initial clusters with KmeansParam. This parameter
  is not used in the final clustering step. iter.max is the maximum
  number of iterations to perform during clustering and sub-clustering
  with KmeansParam. Default value is 30. k is a numeric value indicating
  the k-value used for clustering with NNGraphParam. Default value
  is 10. cluster.fun is a character indicating the graph clustering
  method used with NNGraphParam. By default, the Louvain method is used.

## Value

SpatialExperiment object containing clusters generated from smoothed
data.

## Examples

``` r
data(ClustSignal_example)

# For non-spatial clustering of normalised counts
spe <- clustSIGNAL::p2_clustering(spe)
#> 05:48:54 Calculating PCA using smoothed data.
#> 05:48:54 Final clustering performed on smoothed data. Clusters = 6 
spe$ClustSIGNAL |> head()
#> [1] 3 3 3 3 3 3
#> Levels: 1 2 3 4 5 6
```
