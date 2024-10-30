#' Initial non-spatial clustering
#'
#' @description
#' A function to perform initial non-spatial clustering and sub-clustering of
#' normalised gene expression to generate 'initial clusters'.
#'
#' @param spe SpatialExperiment object. For reclust = FALSE, the object should
#' contain logcounts and PCA, but for reculst = TRUE, the object should contain
#' smoothed gene expression.
#' @param dimRed a character indicating the name of the reduced dimensions to
#' use from the SpatialExperiment object (i.e., from reducedDimNames(spe)).
#' Default value is 'PCA'.
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
#' @return SpatialExperiment object containing 'initial cluster' groups
#' allotted to each cell.
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
#' spe <- clustSIGNAL::p1_clustering(spe, dimRed = "PCA", batch = FALSE,
#'                                  batch_by = "None",
#'                                  clustParams = list(0, 0, 30, 5, "louvain"))
#' head(spe$nsCluster)
#' head(spe$initCluster)
#'
#' @export

#### Non-spatial clustering
p1_clustering <- function(spe, dimRed, batch, batch_by, clustParams) {
    if (clustParams[[1]] == 0) {
        # number of centers = 1/5th of total cells in sample
        clustVal <- min(as.integer(ncol(spe) / 5), 5000)
    } else {
        clustVal <- clustParams[[1]]
    }
    # Initial clustering
    if (batch == TRUE) {
        emb <- harmony::RunHarmony(data_mat = reducedDim(spe, dimRed),
                                   meta_data = colData(spe),
                                   vars_use = batch_by, max.iter = 20,
                                   verbose = FALSE)
        mat <- emb
    } else {
        mat <- reducedDim(spe, dimRed)
    }
    # nsClust <- bluster::clusterRows(
    #     mat,
    #     bluster::TwoStepParam(
    #         first = bluster::KmeansParam(centers = clustVal, iter.max = 30),
    #         second = bluster::NNGraphParam(k = 5, cluster.fun = "louvain")))
    nsClust <- bluster::clusterRows(mat, bluster::TwoStepParam(
        first = bluster::KmeansParam(centers = clustVal,
                                     iter.max = clustParams[[3]]),
        second = bluster::NNGraphParam(k = clustParams[[4]],
                                       cluster.fun = clustParams[[5]])))
    spe$nsCluster <- factor(nsClust)
    show(paste("Initial nonspatial clustering performed. Clusters =",
               length(unique(nsClust)), "Time",
               format(Sys.time(),'%H:%M:%S')))
    # Initial subclustering
    clusters <- length(unique(spe$nsCluster))
    subclusters_list <- list()
    for (cl in seq_len(clusters)){
        speX <- spe[, spe$nsCluster == cl]
        speX <- scater::runPCA(speX)
        matX <- reducedDim(speX, "PCA")
        if (clustParams[[2]] == 0){
            # number of centers = half of total cells in cluster or 1,
            # if only one cell in cluster
            subclustVal <- max(as.integer(ncol(speX) / 2), 1)
        } else {
            subclustVal <- clustParams[[2]]
        }
        # nsClustX <- bluster::clusterRows(
        #     matX,
        #     bluster::TwoStepParam(
        #         first = bluster::KmeansParam(centers = subclustVal,
        #                                      iter.max = 30),
        #         second = bluster::NNGraphParam(k = 5,
        #                                        cluster.fun = "louvain")))
        nsClustX <- bluster::clusterRows(
            matX,
            bluster::TwoStepParam(
                first = bluster::KmeansParam(centers = subclustVal,
                                             iter.max = clustParams[[3]]),
                second = bluster::NNGraphParam(
                    k = clustParams[[4]], cluster.fun = clustParams[[5]])))
        subclusters_list[[cl]] <- setNames(paste0(cl, ".", nsClustX),
                                          colnames(speX))}
    spe$initCluster <- factor(unlist(subclusters_list)[colnames(spe)])
    show(paste("Nonspatial subclustering performed. Subclusters =",
               length(unique(spe$initCluster)), "Time",
               format(Sys.time(),'%H:%M:%S')))
    return (spe)
}
