# Neighbour cell naming

A function to fetch cell IDs.

## Usage

``` r
.cellName(cell, Clust)
```

## Arguments

- cell:

  a vector of neighbourhood cell indices. The cell indices indicate the
  row number of cells in sample metadata.

- Clust:

  a data frame of initial cluster labels of each cell in the sample.

## Value

a data frame of cell IDs of neighbourhood cells.
