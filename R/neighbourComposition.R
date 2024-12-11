#' Cell neighbourhood composition
#'
#' @description
#' A function to calculate the cell neighbourhood composition of 'putative' cell
#' types.
#'
#' @param arr a vector of 'putative cell type' assignments of each cell in the
#' neighbourhood.
#'
#' @return a table of 'putative cell type' proportions in a neighbourhood.
#'
#' @keywords internal

# function to retrieve cluster proportions of each region
.calculateProp <- function(arr) {
    prop <- prop.table(table(arr))
    # check if the proportions for each region sum up to 1
    if (round(sum(prop)) != 1){
        stop("Proportions are incorrect.", "row =", c, "proportion =",
             sum(arr))}
    return(prop)
}
