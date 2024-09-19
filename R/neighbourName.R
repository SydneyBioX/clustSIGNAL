#' Neighbour cell naming
#'
#' @description
#' A function to fetch cell IDs.
#'
#' @param cell a vector of neighbourhood cell indices. The cell indices indicate
#' the row number of cells in sample metadata.
#' @param Clust a data frame of initial cluster assignments of all cells in the
#' sample.
#'
#' @return a data frame of cell IDs of neighbourhood cells.
#'
#' @keywords internal

# function to retrieve cell names of neighbors
.cellName <- function (cell, Clust) {
    return(rownames(Clust)[cell])
}
