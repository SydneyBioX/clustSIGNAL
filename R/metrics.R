#' Clustering metrics
#'
#' @description
#' A function to measure clustering performance using ARI, NMI, and silhouette width, and to provide entropy statistics of each sample in the dataset.
#'
#' @param spe spatialExperiment object with logcounts, PCA, initial non-spatial clustering, entropy, and smoothed gene expression outputs included.
#' @param samples samples a character vector of sample names to which the cells belong. Length of vector must be equal to the number of cells in spatialExperiment object (i.e. the number of rows in colData(spe)).
#' @param celltypes a character vector of cell type annotation of each cell. Length of vector must be equal to the number of cells in spatialExperiment object (i.e. the number of rows in colData(spe)). Set value to NULL, if celltypes information not available.
#' @param cells a character vector of cell IDs of each cell. Length of vector must be equal to the number of cells in spatialExperiment object (i.e. the number of rows in colData(spe)).
#'
#' @return a list containing two items:
#'
#' 1. metrics_out, a table of clustering metrics ARI, NMI, and ASW, and domainness statistics min, max, and mean entropy values calculated for each sample.
#'
#' 2. spe_out, spatialExperiment object with annotated cell type and post-smoothing cluster silhouette widths of each cell added as metadata.
#'
#' @examples

.metrics <- function(spe, samples, celltypes, cells){
    samplesList <- unique(samples)

    if (is.null(celltypes) == TRUE) {
        # Silhoutte width
        silWidthRC <- matrix(nrow = 0, ncol = 3)
        for (s in samplesList) {
            speX <- spe[, samples == s]
            clust_sub <- as.numeric(as.character(speX$reCluster))
            cXg <- t(as.matrix(logcounts(speX)))
            set.seed(12997)
            distMat <- distances(cXg)
            silCluster <- as.matrix(silhouette(clust_sub, distMat))
            silWidthRC <- rbind(silWidthRC, silCluster)
        }
        # QC check
        check.cluster <- identical(silWidthRC[,1], as.numeric(as.character(spe$reCluster)))
        if (check.cluster == FALSE) {
            print(paste("ERROR: issue in silhouette width calculation. Order of cells in silhoutte width matrix did not match cell order in spatial experiment object.", Sys.time()))
        } else {
            spe$rcSil <- silWidthRC[, 3]
        }
        print(paste("No celltype annotations available. Only cluster silhouette widths calculated for each sample.", Sys.time()))
        # other metrics
        tmp_df <- data.frame(Samples = samples,
                             Clusters = spe$reCluster,
                             Entropy = spe$entropy,
                             Silhouette = spe$rcSil)
        metrics_per_sample <- tmp_df %>%
            group_by(as.character(tmp_df$Samples)) %>%
            summarise(ari = NULL,
                      nmi = NULL,
                      asw = mean(Silhouette),
                      minEntropy = min(Entropy),
                      maxEntropy = max(Entropy),
                      meanEntropy = mean(Entropy))
        print(paste("No celltype annotations available. Metrics per sample calculated, except ARI and NMI.", Sys.time()))

    } else {
        # Silhoutte width
        spe$cellNum <- transform(colData(spe), cellNum = as.numeric(interaction(celltypes, drop = TRUE)))$cellNum # integer representatives of cell types
        silWidthRC <- matrix(nrow = 0, ncol = 3)
        silWidthCell <- matrix(nrow = 0, ncol = 3)
        for (s in samplesList) {
            speX <- spe[, samples == s]
            clust_sub <- as.numeric(as.character(speX$reCluster))
            cell_sub <- speX$cellNum
            cXg <- t(as.matrix(logcounts(speX)))
            set.seed(12997)
            distMat <- distances(cXg)
            silCluster <- as.matrix(silhouette(clust_sub, distMat))
            silCell <- as.matrix(silhouette(cell_sub, distMat))
            silWidthRC <- rbind(silWidthRC, silCluster)
            silWidthCell <- rbind(silWidthCell, silCell)
        }
        # QC check
        check.cluster <- identical(silWidthRC[,1], as.numeric(as.character(spe$reCluster)))
        check.cell <- identical(silWidthCell[,1], spe$cellNum)
        if (check.cluster == FALSE | check.cell == FALSE) {
            print(paste("ERROR: issue in silhouette width calculation. Order of cells in silhoutte width matrix did not match cell order in spatial experiment object.", Sys.time()))
        } else {
            spe$rcSil <- silWidthRC[, 3]
            spe$cellSil <- silWidthCell[, 3]
        }
        print(paste("Celltype and cluster silhouette widths calculated for each sample.", Sys.time()))
        # other metrics
        tmp_df <- data.frame(Samples = samples, Celltypes = celltypes,
                             Clusters = spe$reCluster, Entropy = spe$entropy,
                             Silhouette = spe$rcSil)
        metrics_per_sample <- tmp_df %>%
            group_by(as.character(tmp_df$Samples)) %>%
            summarise(ari = ARI(Celltypes, Clusters),
                      nmi = NMI(Celltypes, Clusters),
                      asw = mean(Silhouette),
                      minEntropy = min(Entropy),
                      maxEntropy = max(Entropy),
                      meanEntropy = mean(Entropy))
        print(paste("Metrics per sample calculated.", Sys.time()))
    }
    return(list("metrics_out" = metrics_per_sample,
                "spe_out" = spe))
}
