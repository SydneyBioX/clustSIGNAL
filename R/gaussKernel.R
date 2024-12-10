#' Gaussian distribution weights
#'
#' @description
#' A function to generate weights from a Gaussian distribution.
#'
#' @param ed a numeric value for entropy of the neighbourhood.
#' @param NN an integer for the number of neighbourhood cells the function
#' should consider. The value must be greater than or equal to 1. Default value
#' is 30.
#' @param sd a numeric value for standard deviation of Gaussian distribution.
#' @importFrom stats dnorm
#'
#' @return a matrix where each column contains weights related to specific
#' entropy values.

.gauss_kernel <- function(ed, NN, sd) {
    # Distribution from 0 - entropy, with cells in neighbourhood as cut points

    weight_mat <- dnorm(sapply(ed, function(x) seq(0, x, length.out = NN)),
                     sd = sd)
    # colnames(weight_mat) <- paste0("E", ed)
    return (weight_mat)
}
