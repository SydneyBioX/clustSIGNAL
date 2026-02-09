# Gaussian distribution weights

A function to generate weights from a Gaussian distribution.

## Usage

``` r
.gauss_kernel(ed, NN, sd)
```

## Arguments

- ed:

  a numeric vector of entropy values of all cell neighbourhoods.

- NN:

  an integer for the number of neighbourhood cells including the index
  cell.

- sd:

  a numeric value for standard deviation of Gaussian distribution.

## Value

a matrix where the columns contain weights associated with the entropy
values.
