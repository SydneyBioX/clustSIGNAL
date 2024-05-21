#' Smoothing per neighbourhood
#'
#' @description
#' A function to perform a weighted moving average of gene expression in a neighbourhood.
#'
#' @param mat a gene expression matrix with genes as rows and neighbourhood cells as columns.
#' @param weight a column matrix of weights, where number of rows is equal to number of neighbourhood cells.
#'
#' @return a column matrix of smoothed gene expression.
#'
#' @keywords internal

# function to perform weighted moving average smoothing
.smoothedData <- function(mat, weight) {
    # column matrix of weight
    weight <- as.matrix(weight)
    # return weighted moving average
    return((mat %*% weight) / sum(weight))
}
