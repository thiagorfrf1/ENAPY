namesToAdjacencyKey = """
##
#' @title Names to Adjacency Key
#'
#' @description Convert a vector of stringr, representing names of a square matrix, to an adjacency
#'
#' @details Returns a matrix of 2 rows by choose(length(vector), 2) columns
#'
#' @param vector Vector representing the names of a square matrix
#' @param upper_triangle Not Implemented
#'
#' @export
##
namesToAdjacencyKey <- function(vector, upper_triangle = TRUE) {
  upperTriIndices = triIndices(length(vector)) + 1;
  matrix(vector[upperTriIndices], nrow=2)
}
"""