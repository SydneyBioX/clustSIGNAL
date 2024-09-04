#' Non-spatial clustering
#'
#' @description
#' A function containing two steps used at different times in the clustSIGNAL workflow. An initial non-spatial clustering and sub-clustering step (reclust = FALSE) is used to generate groups of ‘putative cell types’, whereas a later non-spatial clustering step (reclust = TRUE) is used to cluster adaptively smoothed gene expression data.
#'
#' @param spe SpatialExperiment object. For reclust = FALSE, the object should contain logcounts and PCA, but for reculst = TRUE, the object should contain smoothed gene expression.
#' @param dimRed a character indicating the name of the reduced dimensions to use from the SpatialExperiment object (i.e., from reducedDimNames(spe)). Default value is 'PCA'.
#' @param reclust a logical parameter handled within the method.
#' @param ... additional parameters for TwoStepParam clustering methods. Include parameters like k for number of nearest neighbours and cluster.fun for selecting community detection method. Default values k = 5, cluster.fun = "louvain".
#'
#' @return SpatialExperiment object containing 'putative cell type' group allotted to each cell (reclust = FALSE) or clusters generated from smoothed data (reclust = TRUE).
#'
#' @examples
#' data(example)
#'
#' # reclust = FALSE for non-spatial clustering of normalised counts
#' # reclust = TRUE for non-spatial clustering of adaptively smoothed counts
#' spe <- nsClustering(spe, dimRed = "PCA", reclust = FALSE)
#' head(spe$nsCluster)
#' head(spe$nsSubcluster)
#'
#' @export

#### Non-spatial clustering
nsClustering <- function(spe, dimRed = "PCA", reclust, ...) {
    # number of centers = 1/5th of total cells in sample
    clustVal <- min(as.integer(ncol(spe) / 5), 50000)
    if (reclust == FALSE) {
        # Initial clustering
        mat <- reducedDim(spe, dimRed)
        nsClust <- bluster::clusterRows(mat,
                                        bluster::TwoStepParam(
                                            first = bluster::KmeansParam(centers = clustVal, iter.max = 30),
                                            second = bluster::NNGraphParam(k = 5, cluster.fun = "louvain")))
        spe$nsCluster <- factor(nsClust)
        print(paste("Initial nonspatial clustering performed. Clusters =", length(unique(nsClust)), Sys.time()))
        # Initial subclustering
        clusters <- length(unique(spe$nsCluster))
        subclusters_list = list()
        for (c in 1:clusters){
            speX <- spe[, spe$nsCluster == c]
            speX <- scater::runPCA(speX)
            matX <- reducedDim(speX, "PCA")
            # number of centers = half of total cells in cluster or 1,
            # if only one cell in cluster
            subclustVal <- max(as.integer(ncol(speX)/2), 1)
            nsClustX <- bluster::clusterRows(matX,
                                             bluster::TwoStepParam(
                                                 first = bluster::KmeansParam(centers = subclustVal, iter.max = 30),
                                                 second = bluster::NNGraphParam(k = 5, cluster.fun = "louvain")))
            subclusters_list[[c]] <- setNames(paste0(c, ".", nsClustX), colnames(speX))
        }
        spe$nsSubcluster <- factor(unlist(subclusters_list)[colnames(spe)])
        print(paste("Nonspatial subclustering performed. Subclusters =", length(unique(spe$nsSubcluster)), Sys.time()))

    } else if (reclust == TRUE) {
        # Clustering adaptively smoothed data
        spe <- scater::runPCA(spe,
                              assay.type = "smoothed",
                              name = "PCA.smooth")
        mat <- reducedDim(spe, "PCA.smooth")
        reClust <- bluster::clusterRows(mat,
                                        bluster::TwoStepParam(
                                            first = bluster::KmeansParam(centers = clustVal, iter.max = 30),
                                            second = bluster::NNGraphParam(k = 5, cluster.fun = "louvain")))
        spe$reCluster <- factor(reClust)
        print(paste("Nonspatial clustering performed on smoothed data. Clusters =", length(unique(reClust)), Sys.time()))
    }
    return (spe)
}
