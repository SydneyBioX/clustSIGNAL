# Exponential distribution weights

A function to generate weights from an exponential distribution.

## Usage

``` r
.exp_kernel(ed, NN, rate)
```

## Arguments

- ed:

  a numeric vector of entropy values of all cell neighbourhoods.

- NN:

  an integer for the number of neighbourhood cells including the index
  cell.

- rate:

  a numeric value for rate of exponential distribution.

## Value

a matrix where the columns contain weights associated with the entropy
values.
