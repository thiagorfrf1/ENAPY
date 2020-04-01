utils_classes = """

#' Re-class matrix as ena.matrix
#'
#' @param x data.frame, data.table, or matrix to extend
#' @param new.class Additional class to extend the matrix with, default: NULL
#'
#' @return Object of same st
#' @export
as.ena.matrix <- function(x, new.class = NULL) {
  class(x) = c(new.class, "ena.matrix", class(x))
  x
}

#' Re-class matrix as ena.metadata
#'
#' @param x data.frame, data.table, or matrix to extend
#'
#' @return Object of same st
#' @export
as.ena.metadata <- function(x) {
  if(is.factor(x)) {
    x = as.character(x)
  }
  class(x) = c("ena.metadata", "character") # This fails in the $.ena.metadata if is extending character, class(x))
  x
}
as.ena.code <- function(x) {
  if(is.factor(x)) {
    x = as.character(x)
  }
  class(x) = c("ena.code", class(x))
  x
}
as.ena.codes <- function(x) {
  if(is.factor(x)) {
    x = as.character(x)
  }
  class(x) = c("ena.codes", class(x))
  x
}
as.ena.co.occurrence <- function(x) {
  if(is.factor(x)) {
    x = as.character(x)
  }
  class(x) = c("ena.co.occurrence", class(x))
  x
}
as.ena.dimension <- function(x) {
  if(is.factor(x)) {
    x = as.character(x)
  }
  class(x) = c("ena.dimension", class(x))
  x
}
"""