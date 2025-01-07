test_that(
    "Expecting error when wrong column name given for samples", {
        data(ClustSignal_example)
        samples <- "wrongID"

        testthat::expect_error(
            suppressWarnings(clustSIGNAL(spe, samples))
        )
    }
)

test_that(
    "Expecting error when spatial coordinates missing", {
        data(ClustSignal_example)
        samples <- "sample_id"
        colData(spe) <- cbind(colData(spe), spatialCoords(spe))
        spatialCoords(spe) <- NULL

        testthat::expect_error(
            suppressWarnings(clustSIGNAL(spe, samples))
        )
    }
)

