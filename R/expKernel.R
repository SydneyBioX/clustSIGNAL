#' Exponential distribution weights
#'
#' @description
#' A function to generate weights from an exponential distribution.
#'
#' @param ed a numeric value for entropy of the neighbourhood.
#' @param NN an integer for the number of neighbourhood cells the function should consider. The value must be greater than or equal to 1. Default value is 30.
#' @param rate a numeric value for rate of exponential distribution.
#'
#' @return a numeric vector of weights.
#'
#' @examples
#'
#' @keywords internal

# function to retrieve exponential distribution weights for neighbors
.exp_kernel <- function(ed, NN, rate) {
    # distribution from 0 to entropy, with cells in smoothing radius as cut points
    cutpoints <- seq(0, ed, length.out = NN)
    return(dexp(cutpoints, rate = rate))
}
