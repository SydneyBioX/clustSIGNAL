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

