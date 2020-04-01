ena_rotate_by_mean = """
###
#' @title ENA Rotate by mean
#'
#' @description Computes a dimensional reduction from a matrix of points such that
#' the first dimension of the projected space passes through the means of two
#' groups in a the original space. Subsequent dimensions of the projected space
#' are computed using ena.svd
#'
#' @param enaset An \code{\link{ENAset}}
#' @param groups A list containing two logical vectors of length \code{nrow(ENA.set$ena.data$units)},
#' where each vector defines whether a unit is in one of the two groups whose means
#' are used to determine the dimensional reduction
#'
#' @export
#' @return \code{\link{ENARotationSet}}
###
ena.rotate.by.mean <- function(enaset, groups) {
  groups <- list(groups)
  groups <- groups[[1]]
  if (length(groups) < 1) {
    stop("Unable to rotate without 2 groups.")
  }

  if (!is(groups[[1]], "list")) {
    groups <- list(groups);
  }

  data <- as.matrix(enaset$line.weights)
  data <- scale(data, scale = F, center = T);

  col <- NULL
  vals <- NULL

  deflated.data <- data;
  i <- 1;
  weights <- matrix(0, nrow = ncol(deflated.data), ncol = length(groups))

  for (group in 1:length(groups)) {
    col <- group
    vals <- groups[[group]]

    col_one_vals <- deflated.data[vals[[1]], ]
    col_two_vals <- deflated.data[vals[[2]], ]
    col_one_means <- colMeans(as.matrix(col_one_vals))
    col_two_means <- colMeans(as.matrix(col_two_vals))
    col_mean_diff <- col_one_means - col_two_means

    col_mean_diff_sq <- col_mean_diff / sqrt(sum(col_mean_diff ^ 2))

    deflated.data <- deflated.data - (
                      deflated.data %*% col_mean_diff_sq
                    ) %*% t(col_mean_diff_sq)

    weights[, i] <- col_mean_diff_sq
    i <- i + 1;
  }
  defalted_data_svd <- orthogonal_svd(deflated.data, weights);

  colnames(defalted_data_svd) <- c(
    paste0("MR", as.character(1:length(groups))),
    paste0("SVD", as.character((length(groups) + 1):(ncol(defalted_data_svd))))
  )
  rownames(defalted_data_svd) <- colnames(as.matrix(enaset$line.weights))

  rotation_set <- ENARotationSet$new(
    node.positions = NULL,
    rotation = defalted_data_svd,
    codes = enaset$codes
  )
  return(rotation_set)
}

orthogonal_svd <- function(data, weights) {
  if (!is(data, "matrix")) {
    message("orthogonalSVD:  converting data to matrix")
    data <- as.matrix(data)
  }

  #Find the orthogonal transformation that includes W
  Q <- qr_ortho(weights)
  X.bar <- data %*% Q[, (ncol(weights) + 1):ncol(Q)]
  V <- prcomp(X.bar, scale. = F)$rotation
  if (is(V, "numeric")) {
    message("orthogonalSVD:  converting pca result to matrix")
    V <- matrix(V, nrow = length(V))
  }

  to_return <- (cbind(
    Q[, 1:ncol(weights)],
    Q[, (ncol(weights) + 1):ncol(Q)] %*% V
  ))

  return(to_return)
}

qr_ortho <- function(A) {
  return(qr.Q(qr(A), complete = T))
}

"""