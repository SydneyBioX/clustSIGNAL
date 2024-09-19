#' Neighbour cell sorting
#'
#' @description
#' A function to perform neighbourhood cell sorting. Neighbourhood cells that
#' belong to the same 'putative cell type' as the central cell are moved closer
#' to the central cell.
#'
#' @param cell a vector of neighbourhood cell indices. The cell indices indicate
#' the row number of cells in sample metadata.
#' @param Clust a data frame of initial cluster assignments of all cells in the
#' sample.
#'
#' @return a data frame of cell IDs of sorted neighbourhood cells.
#'
#' @keywords internal

# function to sort neighbours and retrieving cell names of neighbors
.cellNameSort <- function (cell, Clust) {
    cellDF <- data.frame(cellName = rownames(Clust)[cell],
                         cluster = Clust[cell, 1])
    # all cells belong to same initial cluster
    if (length(unique(Clust[cell,])) == 1) {
        return(rownames(Clust)[cell])
    } else if (cellDF$cluster[1] %in% cellDF$cluster[seq(2, length(cell))]) {
        centralClust <- cellDF$cluster[1]
        sameClust <- c()
        diffClust <- c()
        for (c in seq_len(nrow(cellDF))) {
            if (cellDF$cluster[c] == centralClust){
                sameClust <- append(sameClust, cellDF$cellName[c])
            } else {
                diffClust <- append(diffClust, cellDF$cellName[c])}}
        region <- c(sameClust, diffClust)
        return(region)
    } else {
        return(rownames(Clust)[cell])}
}
