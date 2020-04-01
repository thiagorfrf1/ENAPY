connection_matrix = """
#' Connection counts as square matrix
#'
#' @param x ena.set or ena.connections (i.e. set$connection.counts)
#'
#' @return matrix
#' @export
connection.matrix <- function(x) {
  if(is(x, "ena.set")) {
    connections <- x$connection.counts
  } else {
    connections <- x
  }
  if(!is(connections, "ena.connections")) {
    stop("Unable to find connections. `x` must be connections from an ena.set or an ena.set")
  }

  simplify <-  (nrow(connections) == 1)
  cm <- as.matrix(connections, square = T, simplify = simplify)
  if(simplify == FALSE && is.list(cm))
    names(cm) <- connections$ENA_UNIT

  return(cm);
}
"""