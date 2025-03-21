#' Cell neighbourhood detection
#'
#' @description
#' A function to identify the neighbourhood of each cell. If sort = TRUE, the
#' neighbourhoods are also sorted such that cells belonging to the same
#' 'initial cluster' as the index cell are arranged closer to it.
#' @param spe SpatialExperiment object with initial cluster and subcluster
#' labels.
#' @param samples a character indicating name of colData(spe) column containing
#' sample names.
#' @param NN an integer for the number of neighbouring cells the function
#' should consider. The value must be greater than or equal to 1. Default value
#' is 30.
#' @param sort a logical parameter for whether to sort the neighbourhood by
#' initial clusters. Default value is TRUE.
#' @param threads a numeric value for the number of CPU cores to be used for the
#' analysis. Default value set to 1.
#' @return a list containing two items:
#'
#' 1. nnCells, a character matrix of NN nearest neighbours - rows are index
#' cells and columns are their nearest neighbours ranging from closest to
#' farthest neighbour. For sort = TRUE, the neighbours belonging to the same
#' initial cluster as the index cell are moved closer to it.
#'
#' 2. regXclust, a list of vectors of each cell's neighbourhood composition
#' indicated by the proportion of initial subclusters it contains.
#' @importFrom BiocNeighbors findKNN KmknnParam
#' @importFrom SpatialExperiment spatialCoords
#' @importFrom methods show
#' @examples
#' data(ClustSignal_example)
#'
#' out_list <- clustSIGNAL::neighbourDetect(spe, samples = "sample_id")
#' out_list |> names()
#' @export

neighbourDetect <- function(spe, samples, NN = 30, sort = TRUE, threads = 1) {
    samplesList <- unique(spe[[samples]])
    nnCells <- matrix(nrow = 0, ncol = NN + 1)
    nnClusts <- matrix(nrow = 0, ncol = NN + 1)
    for (s in samplesList) {
        speX <- spe[, spe[[samples]] == s]
        xy_pos <- spatialCoords(speX)
        Clust <- as.data.frame(as.character(speX$initCluster))
        rownames(Clust) <- colnames(speX)
        subClust <- as.data.frame(speX$initSubcluster)
        rownames(subClust) <- colnames(speX)
        nnMatlist <- BiocNeighbors::findKNN(
            xy_pos, k = NN, num.threads = threads,
            BNPARAM = BiocNeighbors::KmknnParam())
        rownames(nnMatlist$index) <- rownames(subClust)
        nnMatlist$indexNew <- cbind(seq_len(nrow(nnMatlist$index)),
                                    nnMatlist$index)
        if (sort == TRUE) {
            nnCells <- rbind(nnCells, t(apply(nnMatlist$indexNew, 1,
                                              .cellNameSort, Clust = Clust)))
        } else if (sort == FALSE) {
            nnCells <- rbind(nnCells, t(apply(nnMatlist$indexNew, 1,
                                              .cellName, Clust = Clust)))}
        nnClusts <- rbind(nnClusts, t(apply(nnMatlist$indexNew, 1,
                                            .clustNum, subClust = subClust)))}
    regXclust <- apply(nnClusts, 1, .calculateProp)
    # QC check
    check.cells.dim <- identical(dim(nnCells),
                                 as.integer(c(ncol(spe), NN + 1)))
    check.cells.names <- identical(as.character(nnCells[,1]), colnames(spe))
    check.cells.NA <- sum(is.na(nnCells))
    check.clusts.dim <- identical(dim(nnClusts),
                                  as.integer(c(ncol(spe), NN + 1)))
    check.clusts.names <- identical(rownames(nnClusts), colnames(spe))
    check.clusts.NA <- sum(is.na(nnClusts))
    check.region <- identical(length(regXclust), ncol(spe))
    if (check.cells.dim == FALSE) {
        stop("Cell name issue - all neighbors were not detected for all cells.")
    } else if (check.cells.names == FALSE) {
        stop("Cell name issue - Cell order of neighbor data does not match cell
             order in input data.")
    } else if (check.clusts.dim == FALSE) {
        stop("Cluster number issue - all neighbors were not detected for all
             cells.")
    } else if (check.clusts.names == FALSE) {
        stop("Cluster number issues - Cell order of neighbor data does not match
             cell order in input data.")
    } else if (check.region == FALSE) {
        stop("Region proportions not calculated for all cells.")
    } else if (check.cells.NA != 0 | check.clusts.NA != 0) {
        stop("Missing values in neighbor data.")
    } else {
        show(paste("Regions defined. Time", format(Sys.time(),'%H:%M:%S')))}
    return(list("nnCells" = nnCells,
                "regXclust" = regXclust))
}
