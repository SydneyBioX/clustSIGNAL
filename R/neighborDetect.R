#' Cell neighbourhood detection
#'
#' @description
#' A function to identify the neighbourhood of each cell. If sort = TRUE, the neighbourhoods are also sorted such that cells belonging to the same group as the central cell are arranged closer to it.
#'
#' @param spe SpatialExperiment object with logcounts, PCA, and 'putative cell type' groups included.
#' @param samples a character vector of sample names to which the cells belong. Length of vector must be equal to the number of cells in spatialExperiment object (i.e. the number of rows in colData(spe)).
#' @param NN an integer for the number of neighbourhood cells the function should consider. The value must be greater than or equal to 1. Default value is 30.
#' @param cells a character vector of cell IDs of each cell. Length of vector must be equal to the number of cells in SpatialExperiment object (i.e. the number of rows in colData(spe)).
#' @param sort a logical parameter for whether or not to sort the neighbourhood after region description. Default value is TRUE.
#'
#' @return a list containing two items:
#'
#' 1. nnCells, a character matrix of NN nearest neighbours - rows are cells and columns are their nearest neighbours ranged from closest to farthest neighbour. For sort = TRUE, the neighbours belonging to the same 'putative cell type' group as the cell are moved closer to it.
#'
#' 2. regXclust, a list of vectors for each cell's neighbourhood composition indicated by the proportion of 'putative cell type' groups it contains.
#
#' @examples
#'
#' @export

#### Region description + sorting
neighbourDetect <- function(spe, samples, NN, cells, sort) {

    # function to retrieve cell names of neighbors
    cellName <- function (cell, Clust) {
        return(rownames(Clust)[cell])
    }

    # function to sort neighbours and retrieving cell names of neighbors
    cellNameSort <- function (cell, Clust) {
        cellDF = data.frame(cellName = rownames(Clust)[cell], cluster = Clust[cell,1])
        if (cellDF$cluster[1] %in% cellDF$cluster[2:31]) {
            centralClust = cellDF$cluster[1]
            sameClust = c()
            diffClust = c()
            for (c in 1:nrow(cellDF)) {
                if (cellDF$cluster[c] == centralClust){
                    sameClust = append(sameClust, cellDF$cellName[c])
                } else {
                    diffClust = append(diffClust, cellDF$cellName[c])
                }
            }
            region = c(sameClust, diffClust)
            return(region)
        } else {
            return(rownames(Clust)[cell])
        }
    }

    # function to retrieve cluster numbers to which neighbors belong
    clustNum <- function (cell, subClust) {
        return(subClust[cell,])
    }

    # function to retrieve cluster proportions of each region
    calculateProp <- function(arr) {
        prop = prop.table(table(arr))
        if (round(sum(prop)) != 1){ # check if the proportions for each region sum up to 1
            stop(paste("ERROR: proportions are incorrect.", "row =", c, "proportion =", sum(arr)))
        }
        return(prop)
    }

    samplesList <- unique(samples)
    nnCells = matrix(nrow = 0, ncol = NN + 1)
    nnClusts = matrix(nrow = 0, ncol = NN)
    set.seed(12997)
    for (s in samplesList) {
        xy_pos <- spatialCoords(spe[, samples == s])
        Clust <- subset(as.data.frame(colData(spe)), samples == s, nsCluster)
        subClust <- subset(as.data.frame(colData(spe)), samples == s, nsSubcluster)
        nnMatlist <- findKNN(xy_pos, k = NN, BNPARAM = KmknnParam()) # finding NN nearest neighbors for each cell
        rownames(nnMatlist$index) <- rownames(subClust)
        nnMatlist$indexNew <- cbind(seq(1:length(rownames(nnMatlist$index))), nnMatlist$index)
        if (sort == TRUE) {
            nnCells <- rbind(nnCells, t(apply(nnMatlist$indexNew, 1, cellNameSort, Clust = Clust)))
        } else if (sort == FALSE) {
            nnCells <- rbind(nnCells, t(apply(nnMatlist$indexNew, 1, cellName, Clust = Clust)))
        }
        nnClusts <- rbind(nnClusts, t(apply(nnMatlist$index, 1, clustNum, subClust = subClust)))
    }
    # obtaining a list of regions by cluster proportions
    regXclust <- apply(nnClusts, 1, calculateProp)

    # QC check
    check.cells.dim <- identical(dim(nnCells), as.integer(c(length(cells), NN + 1)))
    check.cells.names <- identical(as.character(nnCells[,1]), cells)
    check.cells.NA <- sum(is.na(nnCells))
    check.clusts.dim <- identical(dim(nnClusts), as.integer(c(length(cells), NN)))
    check.clusts.names <- identical(rownames(nnClusts), cells)
    check.clusts.NA <- sum(is.na(nnClusts))
    check.region <- identical(length(regXclust), length(cells))
    if (check.cells.dim == FALSE) {
        stop("ERROR: issue in cell names from neighbour detection. All neighbors were not detected for all cells.")
    } else if (check.cells.names == FALSE) {
        stop("ERROR: issue in cell names from neighbour detection. Order of cells in neighbor data did not match cell order in spatial experiment object.")
    } else if (check.clusts.dim == FALSE) {
        stop("ERROR: issue in cluster numbers from neighbour detection. All neighbors were not detected for all cells.")
    } else if (check.clusts.names == FALSE) {
        stop("ERROR: issue in cluster numbers from neighbour detection. Order of cells in neighbor data did not match cell order in spatial experiment object.")
    } else if (check.region == FALSE) {
        stop("ERROR: Region proportions not calculated for all cells.")
    } else if (check.cells.NA != 0 | check.clusts.NA != 0) {
        stop("ERROR: Missing values in neighbor data.")
    } else {
        print(paste("Regions defined.", Sys.time()))
    }
    return(list("nnCells" = nnCells,
                "regXclust" = regXclust))
}
