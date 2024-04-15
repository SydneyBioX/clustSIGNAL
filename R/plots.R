#' Visualisation plots
#'
#' @description
#' A function to visualise the outputs of clustSIGNAL, including plots for domainness spread and distribution, spatial plots of annotated cells and final clusters, and heatmap comparison of annotated cell types and final clusters.
#'
#' @param spe SpatialExperiment object. For plotType = domain, the object should contain entropy values, for plotType = spatial, the object should contain cluster outputs, and for plotType = heatmap, the object should contain silhouette widths and cluster outputs.
#' @param samples a character vector of sample names to which the cells belong. Length of vector must be equal to the number of cells in SpatialExperiment object (i.e. the number of rows in colData(spe)).
#' @param celltypes a character vector of cell type annotation of each cell. Length of vector must be equal to the number of cells in SpatialExperiment object (i.e. the number of rows in colData(spe)). Set value to NULL, if celltypes information not available.
#' @param plotType a character parameter handled within the method, to decide what plots to generate.
#'
#' @return generates plots depending on the plotType parameter:
#'
#' plotType = domain,  generates a histogram of entropy values and number of cell neighbourhoods and a spatial plot showing entropy of all cell neighbourhoods.
#'
#' plotType = spatial, generates spatial plots of annotated cell types and clusters generated from smoothed data.
#'
#' plotType = heatmap, generates a heat map showing cell number distributions of annotated cell types and post-smoothing clusters.
#'
#' @examples
#'
#' @export

.plots <- function(spe, samples, celltypes, plotType){
    xcoord <- spatialCoords(spe)[, 1]
    ycoord <- spatialCoords(spe)[, 2]
    df <- cbind(as.data.frame(colData(spe)), xcoord, ycoord)
    df$samples <- samples
    midVal <- (max(df$entropy) + min(df$entropy))/2
    entropyRange <- as.integer(max(df$entropy)) + 1

    if (plotType == "domain") {
        # Histogram - entropy spread
        hst <- ggplot(df, aes(entropy)) +
            geom_histogram(binwidth = 0.05) +
            facet_wrap(vars(samples)) +
            ggtitle("Entropy spread of regions") +
            labs(x = "Entropy", y = "Number of regions") +
            theme_grey() +
            theme(text = element_text(size = 20))

        # Spatial coordinates plot - entropy distribution
        sptEntropy <- ggplot(df, aes(x = xcoord, y = -ycoord)) +
            geom_point(size = 1, aes(colour = entropy)) +
            scale_colour_gradient2("Entropy", low = "orangered", mid = "oldlace", high = "blue", midpoint = midVal) +
            scale_size_continuous(range = c(0, entropyRange)) +
            facet_wrap(vars(samples), scales = "free") +
            ggtitle("Spatial distribution of region entropy") +
            labs(x = "x-coordinate", y = "y-coordinate") +
            theme_classic() +
            theme(text = element_text(size = 10))
        print(paste("Domainness plots generated.", Sys.time()))

        return(list("spread" = hst,
                    "distribution" = sptEntropy))


    } else if (plotType == "spatial") {
        sptCluster <- ggplot(df, aes(x = xcoord, y = -ycoord)) +
            geom_point(size = 1, aes(colour = reCluster)) +
            scale_color_manual(values = .getColours()) +
            facet_wrap(vars(samples), scales = "free") +
            ggtitle("Spatial distribution of clusters") +
            labs(x = "x-coordinate", y = "y-coordinate") +
            guides(color = guide_legend(title = "Clusters", override.aes = list(size = 3))) +
            theme_classic() +
            theme(text = element_text(size = 10))

        if (is.null(celltypes) == TRUE) {
            sptCell = NULL
            print(paste("No celltype annotations available. Only cluster spatial plot generated.", Sys.time()))
        } else {
            sptCell <- ggplot(df, aes(x = xcoord, y = -ycoord)) +
                geom_point(size = 1, aes(colour = celltypes)) +
                scale_color_manual(values = .getColours()) +
                facet_wrap(vars(samples), scales = "free") +
                ggtitle("Spatial distribution of annotated cell types") +
                labs(x = "x-coordinate", y = "y-coordinate") +
                guides(color = guide_legend(title = "Cell types", override.aes = list(size = 5))) +
                theme_classic() +
                theme(text = element_text(size = 10))
            print(paste("Celltype and cluster spatial plots generated.", Sys.time()))
        }

        return(list("cellPlot" = sptCell,
                    "clusterPlot" = sptCluster))


    } else if (plotType == "heatmap") {
        # Cluster by cell type heat map
        if (is.null(celltypes) == TRUE) {
            cellHeat = NULL
            print(paste("No celltype annotations available. Celltype and cluster comparison heatmap not generated.", Sys.time()))
        } else {
            meanSilRow <- aggregate(spe$cellSil, list(celltypes), FUN = mean)
            meanSilCol <- aggregate(spe$rcSil, list(spe$reCluster), FUN = mean)
            cellMatrix <- as.matrix(table(celltypes, spe$reCluster))
            colnames(cellMatrix) <- paste("Cluster", colnames(cellMatrix))
            col_func <- colorRamp2(seq(min(cellMatrix), max(cellMatrix), length = 4), c("white", "yellow", "red", "darkred"))
            row_ha <- rowAnnotation(SilhouetteWidth = anno_barplot(meanSilRow[,2], baseline = 0, width = unit(2, "cm")))
            col_ha <- HeatmapAnnotation(SilhouetteWidth = anno_barplot(meanSilCol[,2], baseline = 0, height = unit(2, "cm")))
            cellHeat <- Heatmap(cellMatrix, name = "Number of cells",
                               column_dend_height = unit(4, "cm"),
                               row_dend_width = unit(4, "cm"),
                               top_annotation = col_ha,
                               right_annotation = row_ha,
                               col = col_func,
                               # cluster_rows = FALSE,
                               border_gp = gpar(col = "black", lty = 2),
                               rect_gp = gpar(col = "white", lwd = 2),
                               column_title = "Clusters",
                               column_title_side = "bottom",
                               row_title = "Annotated cell types",
                               cell_fun = function(j, i, x, y, width, height, fill) {
                                   grid.text(sprintf("%.0f", cellMatrix[i, j]), x, y,
                                             gp = gpar(fontsize = 15))
                               })
            print(paste("Celltype and cluster comparison heatmap generated.", Sys.time()))
        }
        return(cellHeat)
    }
}
