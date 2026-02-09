# Adaptive smoothing

A function to perform a weighted, adaptive smoothing of the gene
expression of each cell based on the heterogeneity of its neighbourhood.
Heterogeneous neighbourhoods are smoothed less with higher weights given
to cells belonging to same initial cluster as the index cell.
Homogeneous neighbourhoods are smoothed more with similar weights given
to most cells.

## Usage

``` r
adaptiveSmoothing(spe, nnCells, NN = 30, kernel = "G", spread = 0.3)
```

## Arguments

- spe:

  SpatialExperiment object containing neighbourhood entropy values of
  each cell.

- nnCells:

  a character matrix of NN nearest neighbours - rows are index cells and
  columns are their nearest neighbours ranging from closest to farthest
  neighbour. For sort = TRUE, the neighbours belonging to the same
  initial cluster as the index cell are moved closer to it.

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

## Value

SpatialExperiment object including smoothed gene expression as an
additional assay.

## Examples

``` r
data(ClustSignal_example)

# requires matrix containing NN nearest neighbour cell labels (nnCells),
# generated using the neighbourDetect() function
spe <- clustSIGNAL::adaptiveSmoothing(spe, nnCells)
#> [1] "Smoothing performed. NN = 30 Kernel = G Spread = 0.3 Time 01:12:23"
spe
#> class: SpatialExperiment 
#> dim: 351 1000 
#> metadata(0):
#> assays(2): logcounts smoothed
#> rownames(351): Abcc4 Acp5 ... Zfp57 Zic3
#> rowData names(0):
#> colnames(1000): embryo2_Pos29_cell110_z2 embryo2_Pos29_cell117_z2 ...
#>   embryo2_Pos50_cell361_z5 embryo2_Pos50_cell372_z2
#> colData names(5): uniqueID sample_id entropy initCluster initSubcluster
#> reducedDimNames(1): PCA
#> mainExpName: NULL
#> altExpNames(0):
#> spatialCoords names(2) : X Y
#> imgData names(1): sample_id
```
