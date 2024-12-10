#' Final non-spatial clustering
#'
#' @description
#' A function containing to cluster adaptively smoothed gene expression data.
#'
#' @param spe SpatialExperiment object. For reclust = FALSE, the object should
#' contain logcounts and PCA, but for reculst = TRUE, the object should contain
#' smoothed gene expression.
#' @param batch a logical parameter for whether or not to perform batch
#' correction. Default value is FALSE.
#' @param batch_by a character indicating name of colData(spe) column containing
#' the batch names.
#' @param clustParams a list of parameters for TwoStepParam clustering methods.
#' The clustering parameters are in the order - centers (centers) for clustering
#' with KmeansParam, centers (centers) for sub-clustering clusters with
#' KmeansParam, maximum iterations (iter.max) for clustering with KmeansParam,
#' k values (k) for clustering with NNGraphParam, and community detection method
#' (cluster.fun) to use with NNGraphParam.
#'
#' @return SpatialExperiment object containing clusters generated from smoothed
#' data.
#'
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
#' # For non-spatial clustering of normalised counts
#' spe <- clustSIGNAL::p2_clustering(spe)
#' head(spe$clustSIGNAL)
#'
#' @export

#### Non-spatial clustering
p2_clustering <- function(spe, batch = FALSE, batch_by = "None",
                          clustParams = list(clust_c = 0, subclust_c = 0,
                                             iter.max = 30, k = 5,
                                             cluster.fun = "louvain")) {
    if (clustParams[[1]] == 0){
        # number of centers = 1/5th of total cells in sample
        clustVal <- min(as.integer(ncol(spe) / 5), 5000)
    } else {
        clustVal <- clustParams[[1]]
    }
    # Clustering adaptively smoothed data
    spe <- scater::runPCA(spe, assay.type = "smoothed", name = "PCA.smooth")
    if (batch == TRUE) {
        emb <- harmony::RunHarmony(data_mat = reducedDim(spe, "PCA.smooth"),
                                   meta_data = colData(spe),
                                   vars_use = batch_by, max.iter = 20,
                                   verbose = FALSE)
        mat <- emb
    } else {
        mat <- reducedDim(spe, "PCA.smooth")}
    # reClust <- bluster::clusterRows(mat,
    #                                 bluster::TwoStepParam(
    #                                     first = bluster::KmeansParam(centers = clustVal, iter.max = 30),
    #                                     second = bluster::NNGraphParam(k = 5, cluster.fun = "louvain")))
    reClust <- bluster::clusterRows(
        mat,
        bluster::TwoStepParam(
            first = bluster::KmeansParam(centers = clustVal,
                                         iter.max = clustParams[[3]]),
            second = bluster::NNGraphParam(k = clustParams[[4]],
                                           cluster.fun = clustParams[[5]])))
    spe$clustSIGNAL <- factor(reClust)
    show(paste("Nonspatial clustering performed on smoothed data. Clusters =",
               length(unique(reClust)), "Time",
               format(Sys.time(),'%H:%M:%S')))
    return (spe)
}
