#' Gaussian distribution weights
#'
#' @description
#' A function to generate weights from a Gaussian distribution.
#'
#' @param ed a numeric value for entropy of the neighbourhood.
#' @param NN an integer for the number of neighbourhood cells the function should consider. The value must be greater than or equal to 1. Default value is 30.
#' @param sd a numeric value for standard deviation of Gaussian distribution.
#' @importFrom stats dnorm
#' @return a matrix where each column contains weights related to specific entropy values.

.gauss_kernel <- function(ed, NN, sd) {
    # distribution from 0 to entropy, with cells in smoothing radius as cut points

    # cutpoints <- seq(0, ed, length.out = NN)
    # return(as.matrix(dnorm(cutpoints, sd = sd)))
    weights <- dnorm(sapply(ed, function(x) seq(0, x, length.out = NN)), sd = sd)
    colnames(weights) <- paste0("E", ed)
    return (weights)
}
