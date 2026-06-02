# Heterogeneity measure

A function to measure the heterogeneity of a cell's neighbourhood in
terms of entropy. Generally, homogeneous neighbourhoods have low entropy
and heterogeneous neighbourhoods have high entropy.

## Usage

``` r
entropyMeasure(spe, regXclust)
```

## Arguments

- spe:

  SpatialExperiment object with initial cluster and subcluster labels.

- regXclust:

  a numeric matrix of cells by subclusters, where the values are the
  proportion of initial subclusters in each cell's neighbourhood.

## Value

SpatialExperiment object with entropy values associated with each cell.

## Examples

``` r
data(ClustSignal_example)

# requires matrix containing cluster proportions of each neighbourhood
# (regXclust), generated using the neighbourDetect() function
spe <- clustSIGNAL::entropyMeasure(spe, regXclust)
#> 06:35:23 Neighbourhood heterogeneity calculated.
spe$entropy |> head()
#> embryo2_Pos29_cell110_z2 embryo2_Pos29_cell117_z2 embryo2_Pos29_cell128_z2 
#>                0.4586858                0.2055925                0.6373875 
#> embryo2_Pos29_cell134_z2  embryo2_Pos29_cell14_z2 embryo2_Pos29_cell141_z2 
#>                0.3451173                0.6838104                0.3451173 
```
