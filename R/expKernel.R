#' Exponential distribution weights
#'
#' @description
#' A function to generate weights from an exponential distribution.
#'
#' @param ed a numeric value for entropy of the neighbourhood.
#' @param NN an integer for the number of neighbourhood cells the function
#' should consider. The value must be greater than or equal to 1. Default value
#' is 30.
#' @param rate a numeric value for rate of exponential distribution.
#' @importFrom stats dexp
#'
#' @return a matrix where each column contains weights related to specific
#' entropy values.

.exp_kernel <- function(ed, NN, rate) {
    # distribution from 0 - entropy, with cells in neighbourhood as cut points
    weights <- dexp(sapply(ed, function(x) seq(0, x, length.out = NN)),
                    rate = rate)
    colnames(weights) <- paste0("E", ed)
    return (weights)
}
