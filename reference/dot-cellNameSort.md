# Neighbour cell sorting

A function to perform neighbourhood cell sorting. Neighbourhood cells
that belong to the same initial cluster as the index cell are moved
closer to the index cell.

## Usage

``` r
.cellNameSort(cell, Clust)
```

## Arguments

- cell:

  a vector of neighbourhood cell indices. The cell indices indicate the
  row number of cells in sample metadata.

- Clust:

  a data frame of initial cluster labels of each cell in the sample.

## Value

a data frame of cell IDs of sorted neighbourhood cells.
