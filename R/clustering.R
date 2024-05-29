#' Non-spatial clustering
#'
#' @description
#' A function containing two steps used at different times in the clustSIGNAL workflow. An initial non-spatial clustering and sub-clustering step (reclust = FALSE) is used to generate groups of ‘putative cell types’, whereas a later non-spatial clustering step (reclust = TRUE) is used to cluster adaptively smoothed gene expression data.
#'
#' @param spe SpatialExperiment object. For reclust = FALSE, the object should contain logcounts and PCA, but for reculst = TRUE, the object should contain smoothed gene expression.
#' @param reclust a logical parameter handled within the method.
#'
#' @return SpatialExperiment object containing 'putative cell type' group allotted to each cell (reclust = FALSE) or clusters generated from smoothed data (reclust = TRUE).
#'
#' @examples
#'
#' # to do
#'
#' @export

#### Non-spatial clustering
nsClustering <- function(spe, dimRed, reclust, ...) {
    # number of centers = 1/5th of total cells in sample
    clustVal <- min(as.integer(ncol(spe) / 5), 50000)
    if (reclust == FALSE) {
        if (dimRed == "None") {
            spe <- runPCA(spe)
        }
        else {
            # Clustering
            mat <- reducedDim(spe, dimRed)
            set.seed(12997)
            nsClust <- clusterRows(mat, TwoStepParam(first = KmeansParam(centers = clustVal, iter.max = 30),
                                                     second = NNGraphParam(k = 5, cluster.fun = "louvain")))
            spe$nsCluster <- factor(nsClust)
            print(paste("Initial nonspatial clustering performed. Clusters =", length(unique(nsClust)), Sys.time()))
            # Subclustering
            clusters <- length(unique(spe$nsCluster))
            subclusters_list = list()
            set.seed(12997)
            for (c in 1:clusters){
                speX <- spe[, spe$nsCluster == c]
                speX <- runPCA(speX)
                matX <- reducedDim(speX, "PCA")
                # number of centers = half of total cells in cluster
                subclustVal <- as.integer(ncol(speX) / 2)
                nsClustX <- clusterRows(matX, TwoStepParam(first = KmeansParam(centers = subclustVal, iter.max = 30),
                                                           second = NNGraphParam(k = 5, cluster.fun = "louvain")))
                subclusters_list[[c]] <- setNames(paste0(c, ".", nsClustX), colnames(speX))
            }
            spe$nsSubcluster <- unlist(subclusters_list)[colnames(spe)]
            print(paste("Nonspatial subclustering performed. Subclusters =", length(unique(spe$nsSubcluster)), Sys.time()))
        }
    } else if (reclust == TRUE) {
        # Reclustering
        set.seed(12997)
        spe <- runPCA(spe,
                      assay.type = "smoothed",
                      name = "PCA.smooth")
        mat <- reducedDim(spe, "PCA.smooth")
        set.seed(12997)
        reClust <- clusterRows(mat, TwoStepParam(first = KmeansParam(centers = clustVal, iter.max = 30),
                                                 second = NNGraphParam(k = 5, cluster.fun = "louvain")))
        spe$reCluster <- factor(reClust)
        print(paste("Nonspatial clustering performed on smoothed data. Clusters =", length(unique(reClust)), Sys.time()))
    }
    return (spe)
}
