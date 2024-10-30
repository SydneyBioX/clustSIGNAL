# test_that("Checking dimensions of neighbours matrix", {
#     data(example)
#     cells <- "uniqueID"
#     samples <- "sample_id"
#     NN <- 30
#     out <- clustSIGNAL::neighbourDetect(spe, samples, NN, cells,
#                                         sort = TRUE)
#
#     testthat::expect_equal(
#         dim(out$nnCells),
#         as.integer(c(ncol(spe), NN + 1))
#     )
# })
#
# test_that("Checking neighbours detected for all cells", {
#     data(example)
#     cells <- "uniqueID"
#     samples <- "sample_id"
#     NN <- 30
#     out <- clustSIGNAL::neighbourDetect(spe, samples, NN, cells,
#                                         sort = TRUE)
#
#     testthat::expect_identical(
#         as.character(out$nnCells[, 1]),
#         spe[[cells]]
#     )
# })

test_that(
    "Expecting error when wrong column name given for samples", {
        data(example)
        samples <- "wrongID"
        cells <- "uniqueID"

        testthat::expect_error(
            suppressWarnings(clustSIGNAL(spe, samples, cells))
        )
    }
)

test_that(
    "Expecting error when wrong column name given for cells", {
        data(example)
        samples <- "sample_id"
        cells <- "wrongID"

        testthat::expect_error(
            suppressWarnings(clustSIGNAL(spe, samples, cells))
        )
    }
)

test_that(
    "Expecting error when spatial coordinates missing", {
        data(example)
        samples <- "sample_id"
        cells <- "uniqueID"
        colData(spe) <- cbind(colData(spe), spatialCoords(spe))
        spatialCoords(spe) <- NULL

        testthat::expect_error(
            suppressWarnings(clustSIGNAL(spe, samples, cells))
        )
    }
)

