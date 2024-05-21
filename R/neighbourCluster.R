#' Neighbour cell 'putative cell type' assignment
#'
#' @description
#' A function to fetch 'putative cell type' assignments of cell.
#'
#' @param cell a vector of neighbourhood cell indices. The cell indices indicate the row number of cells in sample metadata.
#' @param subClust a data frame of 'putative cell type' assignments of all cells in the sample.
#'
#' @return a data frame of 'putative cell type' assignments of neighbourhood cells.
#'
#' @keywords internal

# function to retrieve cluster numbers to which neighbors belong
.clustNum <- function (cell, subClust) {
    return(subClust[cell,])
}
