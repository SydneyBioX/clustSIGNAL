# Heterogeneity measure

A function to measure the heterogeneity of a cell's neighbourhood in
terms of entropy. Generally, homogeneous neighbourhoods have low entropy
and heterogeneous neighbourhoods have high entropy.

## Usage

``` r
entropyMeasure(spe, regXclust, threads = 1)
```

## Arguments

- spe:

  SpatialExperiment object with initial cluster and subcluster labels.

- regXclust:

  a list of vectors of each cell's neighbourhood composition indicated
  by the proportion of initial subclusters it contains.

- threads:

  a numeric value for the number of CPU cores to be used for the
  analysis. Default value set to 1.

## Value

SpatialExperiment object with entropy values associated with each cell.

## Examples

``` r
data(ClustSignal_example)

# requires list containing cluster proportions of each region (regXclust),
# generated using the neighbourDetect() function
spe <- clustSIGNAL::entropyMeasure(spe, regXclust)
#> [1] "Region heterogeneity calculated. Time 01:12:27"
spe$entropy |> head()
#> [1] 0 0 0 0 0 0
```
