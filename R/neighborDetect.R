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
#' @importFrom BiocNeighbors findKNN KmknnParam
#' @importFrom SpatialExperiment spatialCoords
#' @importFrom methods show
#'
#' @examples
#' data(example)
#'
#' out_list <- neighbourDetect(spe, samples = "sample_id", NN = 30, cells = "uniqueID", sort = TRUE)
#' names(out_list)
#'
#' @export

#### Region description + sorting
neighbourDetect <- function(spe, samples, NN = 30, cells, sort = TRUE) {
    samplesList <- unique(spe[[samples]])
    nnCells <- matrix(nrow = 0, ncol = NN + 1)
    nnClusts <- matrix(nrow = 0, ncol = NN)
    for (s in samplesList) {
        speX <- spe[, spe[[samples]] == s]
        xy_pos <- spatialCoords(speX)
        Clust <- as.data.frame(as.character(speX$nsCluster))
        rownames(Clust) <- speX[[cells]]
        subClust <- as.data.frame(speX$nsSubcluster)
        rownames(subClust) <- speX[[cells]]
        # finding NN nearest neighbors for each cell
        nnMatlist <- BiocNeighbors::findKNN(xy_pos, k = NN, BNPARAM = BiocNeighbors::KmknnParam())
        rownames(nnMatlist$index) <- rownames(subClust)
        # adding central cell indices as first column of matrix of NN neighborhood indices
        nnMatlist$indexNew <- cbind(seq_len(nrow(nnMatlist$index)), nnMatlist$index)
        if (sort == TRUE) {
            nnCells <- rbind(nnCells, t(apply(nnMatlist$indexNew, 1, .cellNameSort, Clust = Clust)))
        } else if (sort == FALSE) {
            nnCells <- rbind(nnCells, t(apply(nnMatlist$indexNew, 1, .cellName, Clust = Clust)))
        }
        nnClusts <- rbind(nnClusts, t(apply(nnMatlist$index, 1, .clustNum, subClust = subClust)))
    }
    # obtaining a list of regions by cluster proportions
    regXclust <- apply(nnClusts, 1, .calculateProp)

    # QC check
    check.cells.dim <- identical(dim(nnCells), as.integer(c(length(spe[[cells]]), NN + 1)))
    check.cells.names <- identical(as.character(nnCells[,1]), spe[[cells]])
    check.cells.NA <- sum(is.na(nnCells))
    check.clusts.dim <- identical(dim(nnClusts), as.integer(c(length(spe[[cells]]), NN)))
    check.clusts.names <- identical(rownames(nnClusts), spe[[cells]])
    check.clusts.NA <- sum(is.na(nnClusts))
    check.region <- identical(length(regXclust), length(spe[[cells]]))
    if (check.cells.dim == FALSE) {
        stop("Issue in cell names from neighbour detection. All neighbors were not detected for all cells.")
    } else if (check.cells.names == FALSE) {
        stop("Issue in cell names from neighbour detection. Order of cells in neighbor data did not match cell order in spatial experiment object.")
    } else if (check.clusts.dim == FALSE) {
        stop("Issue in cluster numbers from neighbour detection. All neighbors were not detected for all cells.")
    } else if (check.clusts.names == FALSE) {
        stop("Issue in cluster numbers from neighbour detection. Order of cells in neighbor data did not match cell order in spatial experiment object.")
    } else if (check.region == FALSE) {
        stop("Region proportions not calculated for all cells.")
    } else if (check.cells.NA != 0 | check.clusts.NA != 0) {
        stop("Missing values in neighbor data.")
    } else {
        show(paste("Regions defined.", Sys.time()))
    }
    return(list("nnCells" = nnCells,
                "regXclust" = regXclust))
}
