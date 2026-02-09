# Cell neighbourhood detection

A function to identify the neighbourhood of each cell. If sort = TRUE,
the neighbourhoods are also sorted such that cells belonging to the same
'initial cluster' as the index cell are arranged closer to it.

## Usage

``` r
neighbourDetect(spe, samples, NN = 30, sort = TRUE, threads = 1)
```

## Arguments

- spe:

  SpatialExperiment object with initial cluster and subcluster labels.

- samples:

  a character indicating name of colData(spe) column containing sample
  names.

- NN:

  an integer for the number of neighbouring cells the function should
  consider. The value must be greater than or equal to 1. Default value
  is 30.

- sort:

  a logical parameter for whether to sort the neighbourhood by initial
  clusters. Default value is TRUE.

- threads:

  a numeric value for the number of CPU cores to be used for the
  analysis. Default value set to 1.

## Value

a list containing two items:

1\. nnCells, a character matrix of NN nearest neighbours - rows are
index cells and columns are their nearest neighbours ranging from
closest to farthest neighbour. For sort = TRUE, the neighbours belonging
to the same initial cluster as the index cell are moved closer to it.

2\. regXclust, a list of vectors of each cell's neighbourhood
composition indicated by the proportion of initial subclusters it
contains.

## Examples

``` r
data(ClustSignal_example)

out_list <- clustSIGNAL::neighbourDetect(spe, samples = "sample_id")
#> [1] "Regions defined. Time 01:12:28"
out_list |> names()
#> [1] "nnCells"   "regXclust"
```
