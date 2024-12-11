## clustSIGNAL v0.99.6 (2024-12-11)
* Modified adaptiveSmoothing() function to improve its runtime.
* Improved the tutorial vignette content.
* Added step-by-step ClustSIGNAL run guide to the tutorial vignette.
* Removed the user parameter cells. The package automatically uses the column names of the input SpatialExperiment object as the cell ID.
* Added a cell id check - throws error if duplicates found among column names of input SpatialExperiment object.

## clustSIGNAL v0.99.5 (2024-11-28)
* Fixed issue where absence of spatial coordinates was not throwing error at the beginning of the run.

## clustSIGNAL v0.99.4 (2024-11-28)
* Corrected character check for output format type.

## clustSIGNAL v0.00.3 (2024-10-31)
* Clustering parameter options added to package
* Initial and final clustering steps are now in separate functions
* Test units added to the package

## clustSIGNAL v0.00.2 (2024-10-08)
* Updated parallel runs in smoothing function for faster runs

## clustSIGNAL v0.99.0 (2024-09-19)
* Submitted to Bioconductor.
