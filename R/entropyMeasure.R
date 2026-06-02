#' Heterogeneity measure
#'
#' @description
#' A function to measure the heterogeneity of a cell's neighbourhood in terms of
#' entropy. Generally, homogeneous neighbourhoods have low entropy and
#' heterogeneous neighbourhoods have high entropy.
#' @param spe SpatialExperiment object with initial cluster and subcluster
#' labels.
#' @param regXclust a numeric matrix of cells by subclusters, where the values
#' are the proportion of initial subclusters in each cell's neighbourhood.
#' @return SpatialExperiment object with entropy values associated with each
#' cell.
#' @examples
#' data(ClustSignal_example)
#'
#' # requires matrix containing cluster proportions of each neighbourhood
#' # (regXclust), generated using the neighbourDetect() function
#' spe <- clustSIGNAL::entropyMeasure(spe, regXclust)
#' spe$entropy |> head()
#' @export

entropyMeasure <- function(spe, regXclust) {
    # calculating max entropy
    max_ent <- log2(ncol(regXclust))
    nz_prop <- regXclust > 0
    p_log <- matrix(0, nrow = nrow(regXclust), ncol = ncol(regXclust))
    p_log[nz_prop] <- regXclust[nz_prop] * log2(regXclust[nz_prop])
    regEntropy <- setNames(-rowSums(p_log), rownames(regXclust))
    # # normalising entropy by total subclusters
    # if (max_ent > 0) { # i.e., only 1 cluster
    #     norm_entropy <- regEntropy / max_ent
    # } else {
    #     norm_entropy <- regEntropy
    # }
    # norm_entropy <- round(norm_entropy, 5)

    # QC checks
    if (!identical(names(regEntropy), colnames(spe))) {
        stop("Cells have likely been allocated wrong entropy values.")
    }
    if (anyNA(regEntropy)) {
        stop("Entropy values are missing for some cells.")
    }

    spe$entropy <- regEntropy
    message(sprintf("%s Neighbourhood heterogeneity calculated.",
                    format(Sys.time(),'%H:%M:%S')))
    return(spe)
}
