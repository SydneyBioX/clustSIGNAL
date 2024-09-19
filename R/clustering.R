#' Non-spatial clustering
#'
#' @description
#' A function containing two steps used at different times in the clustSIGNAL
#' workflow. An initial non-spatial clustering and sub-clustering step (reclust =
#' FALSE) is used to generate groups of ‘putative cell types’, whereas a later
#' non-spatial clustering step (reclust = TRUE) is used to cluster adaptively
#' smoothed gene expression data.
#'
#' @param spe SpatialExperiment object. For reclust = FALSE, the object should
#' contain logcounts and PCA, but for reculst = TRUE, the object should contain
#' smoothed gene expression.
#' @param samples a character indicating name of colData(spe) column containing
#' sample names.
#' @param dimRed a character indicating the name of the reduced dimensions to
#' use from the SpatialExperiment object (i.e., from reducedDimNames(spe)).
#' Default value is 'PCA'.
#' @param batch a logical parameter for whether or not to perform batch
#' correction. Default value is FALSE.
#' @param reclust a logical parameter handled within the method.
#' @param ... additional parameters for TwoStepParam clustering methods. Include
#' parameters like k for number of nearest neighbours and cluster.fun for
#' selecting community detection method. Default values k = 5, cluster.fun =
#' "louvain".
#'
#' @return SpatialExperiment object containing 'putative cell type' group
#' allotted to each cell (reclust = FALSE) or clusters generated from smoothed
#' data (reclust = TRUE).
#' @importFrom bluster clusterRows TwoStepParam KmeansParam NNGraphParam
#' @importFrom scater runPCA
#' @importFrom harmony RunHarmony
#' @importFrom SummarizedExperiment colData
#' @importFrom SingleCellExperiment reducedDim
#' @importFrom methods show
#' @importFrom stats setNames
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
nsClustering <- function(spe, samples, dimRed = "PCA", batch = FALSE,
                         reclust, ...) {
    # number of centers = 1/5th of total cells in sample
    clustVal <- min(as.integer(ncol(spe) / 5), 5000)
    if (reclust == FALSE) {
        # Initial clustering
        if (batch == TRUE) {
            emb <- harmony::RunHarmony(data_mat = reducedDim(spe, dimRed),
                                       meta_data = colData(spe),
                                       vars_use = samples, max.iter = 20,
                                       verbose = FALSE)
            mat <- emb
        } else {
            mat <- reducedDim(spe, dimRed)}
        nsClust <- bluster::clusterRows(mat, bluster::TwoStepParam(
            first = bluster::KmeansParam(centers = clustVal, iter.max = 30),
            second = bluster::NNGraphParam(k = 5, cluster.fun = "louvain")))
        spe$nsCluster <- factor(nsClust)
        show(paste("Initial nonspatial clustering performed. Clusters =",
                   length(unique(nsClust)), "Time",
                   format(Sys.time(),'%H:%M:%S')))
        # Initial subclustering
        clusters <- length(unique(spe$nsCluster))
        subclusters_list <- list()
        for (c in seq_len(clusters)){
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
            subclusters_list[[c]] <- setNames(paste0(c, ".", nsClustX),
                                              colnames(speX))}
        spe$nsSubcluster <- factor(unlist(subclusters_list)[colnames(spe)])
        show(paste("Nonspatial subclustering performed. Subclusters =",
                   length(unique(spe$nsSubcluster)), "Time",
                   format(Sys.time(),'%H:%M:%S')))
    } else if (reclust == TRUE) {
        # Clustering adaptively smoothed data
        spe <- scater::runPCA(spe, assay.type = "smoothed", name = "PCA.smooth")
        if (batch == TRUE) {
            emb <- harmony::RunHarmony(data_mat = reducedDim(spe, "PCA.smooth"),
                                       meta_data = colData(spe),
                                       vars_use = samples, max.iter = 20,
                                       verbose = FALSE)
            mat <- emb
        } else {
            mat <- reducedDim(spe, "PCA.smooth")}
        reClust <- bluster::clusterRows(mat,
                                        bluster::TwoStepParam(
                                            first = bluster::KmeansParam(centers = clustVal, iter.max = 30),
                                            second = bluster::NNGraphParam(k = 5, cluster.fun = "louvain")))
        spe$clustSIGNAL <- factor(reClust)
        show(paste("Nonspatial clustering performed on smoothed data. Clusters =",
                   length(unique(reClust)), "Time",
                   format(Sys.time(),'%H:%M:%S')))}
    return (spe)
}
