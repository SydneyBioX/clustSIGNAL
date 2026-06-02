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
#' 2. regXclust, a numeric matrix of each cell's neighbourhood composition
#' indicated by the proportion of initial subclusters (column) in each cell
#' (row).
#'
#' @importFrom BiocNeighbors findKNN KmknnParam
#' @importFrom SpatialExperiment spatialCoords
#' @importFrom methods show
#' @importFrom SummarizedExperiment colData
#' @examples
#' data(ClustSignal_example)
#'
#' out_list <- clustSIGNAL::neighbourDetect(spe, samples = "sample_id")
#' out_list |> names()
#' @export

neighbourDetect <- function(spe, samples, NN = 30, sort = TRUE, threads = 1) {
    # extracting sample names
    samplesList <- unique(colData(spe)[[samples]])

    # finding neighbours for each sample
    results_list <- lapply(samplesList, function(s) {
        speX <- spe[, colData(spe)[[samples]] == s]
        xy_pos <- spatialCoords(speX)
        Cell_names <- colnames(speX)
        Clust <- as.character(speX$initCluster)
        subClust <- as.character(speX$initSubcluster)
        nnMatlist <- findKNN(xy_pos, k = NN, num.threads = threads,
                             BNPARAM = KmknnParam())
        rownames(nnMatlist$index) <- colnames(speX)
        # adding index cell number as first column
        indexNew <- cbind(seq_len(nrow(nnMatlist$index)),
                          nnMatlist$index)
        # sorting neighbourhood
        mat_clusts <- matrix(subClust[indexNew], nrow = nrow(indexNew),
                             ncol = ncol(indexNew))
        rownames(mat_clusts) <- colnames(speX)
        if (isTRUE(sort)) {
            # Extract the unsorted names and initial clusters natively
            mat_cells_unsort <- matrix(Cell_names[indexNew],
                                       nrow = nrow(indexNew))
            mat_clusts_init <- matrix(Clust[indexNew], nrow = nrow(indexNew))
            # sorting neighbourhood
            mat_cells <- t(sapply(seq_len(nrow(indexNew)), function(i) {
                names_row  <- mat_cells_unsort[i, ]
                clusts_row <- mat_clusts_init[i, ]
                # find neighbouring cells from same cluster as index cell
                is_match <- clusts_row == clusts_row[1]
                # rearrange neighbours
                return(names_row[c(which(is_match), which(!is_match))])
            }))
        } else {
            mat_cells <- matrix(Cell_names[indexNew], nrow = nrow(indexNew),
                                ncol = ncol(indexNew))
        }
        rownames(mat_cells) <- colnames(speX)
        return(list(cells = mat_cells, clusts = mat_clusts))
    })

    # combining data from all samples
    nnCells  <- do.call(rbind, lapply(results_list, `[[`, "cells"))
    nnClusts <- do.call(rbind, lapply(results_list, `[[`, "clusts"))

    # generating cluster proportions
    all_clusters <- sort(unique(as.vector(nnClusts)))
    nnClusts_int <- matrix(match(nnClusts, all_clusters),
                           nrow = nrow(nnClusts),
                           ncol = ncol(nnClusts))
    regXclust <- t(apply(nnClusts_int, 1, function(row_arr) {
        prop <- tabulate(row_arr, nbins = length(all_clusters)) / length(row_arr)
        # QC check
        if (round(sum(prop)) != 1) {
            stop("Neighbourhood cluster proportions are incorrect.")
        }
        return(prop)
    }))
    rownames(regXclust) <- rownames(nnClusts)
    colnames(regXclust) <- all_clusters

    # QC checks
    target_dim <- as.integer(c(ncol(spe), NN + 1))

    if (!identical(dim(nnCells), target_dim) |
        !identical(dim(nnClusts), target_dim)) {
        stop("Cell name or cluster number issue: neighbors were not detected for
             some cells.")
    }
    if (!identical(as.character(nnCells[, 1]), colnames(spe)) |
        !identical(rownames(nnClusts), colnames(spe))) {
        stop("Cell name or cluster number issue: cell ids are misaligned between
             neighbour data and original input data.")
    }
    if (!identical(nrow(regXclust), ncol(spe))) {
        stop("Neighbourhood composition not assessed for some cells.")
    }
    if (anyNA(nnCells) | anyNA(nnClusts)) {
        stop("Missing values in neighbor data.")
    }

    message(sprintf("%s Neighbourhoods defined.", format(Sys.time(), '%H:%M:%S')))

    return(list("nnCells" = nnCells, "regXclust" = regXclust))
}
