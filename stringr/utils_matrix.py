utils_matrix = """
#' ENA line weights as matrix
#'
#' @param x ena.line.weights data.table to covert to matrix
#' @param ... additional arguments to be passed to or from methods
#' @param square [TBD]
#'
#' @return matrix
#' @export
as.matrix.ena.line.weights <- function(x, ..., square = FALSE) {
  args = list(...)

  # if(!is.null(args$square))
  #   square = args$square

  class(x) = class(x)[-1]
  rows = x[, !find_meta_cols(x), with = F]

  if(square) {
    upperTriSize = ncol(rows)
    number = ( (ceiling(sqrt(2*upperTriSize)) ^ 2) ) - (2*upperTriSize)
    codes = unique(unlist(sapply(colnames(rows), strsplit, split = " & ")))
    cm = sapply(seq(nrow(rows)), function(unit) {
      m = matrix(NA, number,  number, dimnames = list(codes, codes))
      m[upper.tri(m)] = as.numeric(rows[unit,])
      m
    }, simplify = F);
    return(cm)
  }
  else {
   as.matrix(remove_meta_data(rows), ...)
  }
}

#' ENA rotations as matrix
#'
#' @param x ena.rotation.matrix to conver to matrix
#' @param ... 	additional arguments to be passed to or from methods
#'
#' @return matrix
#' @export
as.matrix.ena.rotation.matrix <- function(x, ...) {
  class(x) = class(x)[-1]
  x = remove_meta_data(x)
  as.matrix(x, ...)
}

#' ENA points as matrix
#'
#' @param x ena.points to convert to a matrix
#' @param ... 	additional arguments to be passed to or from methods
#'
#' @return matrix
#' @export
as.matrix.ena.points <- function(x, ...) {
  class(x) = class(x)[-1]
  x = remove_meta_data(x)
  as.matrix(x, ...)
}

#' Matrix without metadata
#'
#' @param x Object to convert to  a matrix
#' @param ... 	additional arguments to be passed to or from methods
#'
#' @return matrix
#' @export
as.matrix.ena.matrix <- function(x, ...) {
  class(x) = class(x)[-1]
  x = remove_meta_data(x)
  as.matrix(x, ...)
}

#' ENA nodes as matrix
#'
#' @param x ena.nodes to convert to matrix
#' @param ... 	additional arguments to be passed to or from methods
#'
#' @return matrix
#' @export
as.matrix.ena.nodes <- function(x, ...) {
  class(x) = class(x)[-1]
  as.matrix(x[,-c("code")], ...)
}

#' ENA row connections as matrix
#'
#' @param x ena.row.connections to conver to a matrix
#' @param ... 	additional arguments to be passed to or from methods
#'
#' @return matrix
#' @export
as.matrix.row.connections <- function(x, ...) {
  class(x) = class(x)[-1]
  as.matrix(x[, sapply(x, is, class2="ena.co.occurrence"), with = F], ...)
}


#' ENA Connections as a matrix
#'
#' @param x ena.connections object
#' @param ... 	additional arguments to be passed to or from methods
# @param square Logical. If TRUE, each row is converted to a square matrix
# @param simplify Logical. If TRUE, returns back a single result as vector
# @param names Ignored
#'
#' @return If square is FALSE (default), a matrix with all metadata columns removed, otherwise a list with square matrices
#' @export
as.matrix.ena.connections <- function(x, ...) {
  class(x) = class(x)[-1]
  xx = remove_meta_data(x)
  rows = as.data.frame(x)[, !find_meta_cols(x), drop = F]

  args = list(...)
  if(is.null(args$square))
    square = F
  else
    square = args$square

  names = args$names

  if(is.null(args$simplify))
    simplify = ifelse(nrow(x) > 1, F, T)
  else
    simplify = args$simplify

  if(square) {
    upperTriSize = ncol(rows)
    number = ( (ceiling(sqrt(2*upperTriSize)) ^ 2) ) - (2*upperTriSize)
    codes = unique(unlist(sapply(colnames(rows), strsplit, split = " & ")))
    cm = sapply(seq(nrow(rows)), function(unit) {
      m = matrix(NA, number, number, dimnames = list(codes, codes))
      m[upper.tri(m)] = as.numeric(rows[unit,])
      m
    }, simplify = F)

    if(simplify) {
      cm = cm[[1]]
    } else {
      names(cm) = names
    }
  } else {
    cm = as.matrix(rows)
    rownames(cm) = names
  }

  cm
}
"""