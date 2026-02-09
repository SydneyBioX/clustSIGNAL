# Initial non-spatial clustering

A function to perform initial non-spatial clustering and sub-clustering
of normalised gene expression to generate 'initial clusters' and
'initial subclusters'.

## Usage

``` r
p1_clustering(
  spe,
  dimRed = "None",
  batch = FALSE,
  batch_by = "None",
  threads = 1,
  clustParams = list(clust_c = 0, subclust_c = 0, iter.max = 30, k = 10, cluster.fun =
    "louvain")
)
```

## Arguments

- spe:

  a SpatialExperiment object containing spatial coordinates in
  'spatialCoords' matrix and normalised gene expression in 'logcounts'
  assay.

- dimRed:

  a character indicating the name of the reduced dimensions to use from
  the SpatialExperiment object (i.e., from reducedDimNames(spe)).
  Default value is 'None'.

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
  default set to 0, in which case the method uses either 5000 centers or
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

SpatialExperiment object with initial cluster and subcluster labels of
each cell.

## Examples

``` r
data(ClustSignal_example)

spe <- clustSIGNAL::p1_clustering(spe, dimRed = "PCA")
#> [1] "Initial nonspatial clustering performed. Clusters = 4 Time 01:12:29"
#> [1] "Nonspatial subclustering performed. Subclusters = 12 Time 01:12:29"
spe$nsCluster |> head()
#> NULL
spe$initCluster |> head()
#> [1] 2 2 3 2 2 2
#> Levels: 1 2 3 4
```
