utils ="""
#' Find metadata columns
#'
#' @param x data.table (or frame) to search for columns of class ena.metadata
#'
#' @return logical vector
#' @export
find_meta_cols <- function(x) {
   sapply(x, is, class2 = "ena.metadata")
}

#' Find code columns
#'
#' @param x data.table (or frame) to search for columns of class ena.co.occurrence
#'
#' @return logical vector
#' @export
find_code_cols <- function(x) {
   grepl("adjacency.code", x = names(x)) | sapply(x, function(col) {
     is(col, class2 = "ena.co.occurrence")
   })
}

#' Find dimension columns
#'
#' @param x data.table (or frame) to search for columns of class ena.dimension
#'
#' @return logical vector
#' @export
find_dimension_cols <- function(x) {
   sapply(x, is, class2 = "ena.dimension")
}

#' Remove meta columns from data.table
#'
#' @param x [TBD]
#'
#' @return data.table withe columns of class ena.meta.data removed
#' @export
remove_meta_data <- function(x) {
   as.data.frame(x)[, !find_meta_cols(x), drop = F]
}

#' Extract metadata easily
#'
#' @param x [TBD]
#' @param i [TBD]
#'
#' @return [TBD]
#' @export
"$.ena.metadata" <- function(x, i) {
   #browser()
   parts <- unlist(strsplit(
               x = as.character(sys.call())[2], split = "\\$"
            ))[1:2]

   set <- get(parts[1], envir = parent.frame())
   m <- set[[parts[2]]][x == i, ]
   m
}

#' Extract line.weignts easily
#'
#' @param x [TBD]
#' @param i [TBD]
#'
#' @return [TBD]
#' @export
"$.line.weights" <- function (x, i) {
   vals <- x[[which(colnames(x) == i)]]

   vals
}

#' Extract points easily
#'
#' @param x [TBD]
#' @param i [TBD]
#'
#' @return [TBD]
#' @export
"$.ena.points" <- function (x, i) {
   vals <- x[[which(colnames(x) == i)]]

   vals
}
# "$.ena.plot" <- function(x, i) {
#  browser()
# }
# "[[.ena.plot" <- function(x, i) {
#  browser()
# }
#' @export
.DollarNames.ena.metadata <- function(x, pattern = "") {
   unique(x)
}

# "[.ena.matrix" = function(x, ...)
# {
#    browser()
#    original.class = class(x)[1]
#    class(x) = class(x)[-1]
#    x = x[...]
#
# #   y = as.data.frame(x)
# }

#' @export
summary.ena.set <- function(object, ...) {
   x <- object
   print_dims <- function(n = 2) {
      cat("\t", paste("Dimension", 1:n, collapse = "\t"), "\n")
   }
   cat("Units: ", nrow(x$points), "\t\t")
   cat("Codes: ", length(x$rotation$codes), "\n")

   cat("Variance: \n")
   print_dims()
   cat("\t", paste(round(x$model$variance[1:2], 3), collapse = "\t\t"), "\n\n")

   cat("Eigenvalues: \n")
   print_dims()
   cat("\t", paste(round(
      x$rotation$eigenvalues[1:2], 3), collapse = "\t\t"), "\n\n")

   cat("Correlations: \n")
   cors <- ena.correlations(x)
   rownames(cors) <- paste("Dimension", 1:2)
   print(cors)
}
# as.data.frame.ena.connections <- function(x) {
#   class(x) = class(x)[-1]
#   y = as.data.frame(x)
#   y
# }
# format.co.occurrence = format.metadata = function(x, justify = "none") {
#   y = as.character(x)
#   format(y, justify = justify)
# }

#' Title
#'
#' @param x [TBD]
#' @param ... [TBD]
#' @param plot [TBD]
#' @param set [TBD]
#'
#' @return [TBD]
#' @export
print.ena.set <- function(x, ..., plot = FALSE, set = TRUE) {
   x.unclass <- unclass(x)

   if(
      !is.null(x.unclass$`_plot_op`) &&
      x.unclass$`_plot_op` == T
   ) {
      base::print(x.unclass$plots)
   }
   else {
      if(plot == FALSE) {
         x.unclass$plots <- NULL
      }
      base::print(x.unclass)
   }
}

#' Title
#'
#' @param x [TBD]
#' @param by [TBD]
#' @param model [TBD]
#' @param ... [TBD]
#'
#' @return [TBD]
#' @export
as_trajectory <- function(x,
   by = x$`_function.params`$conversation[1],
   model = c("AccumulatedTrajectory", "SeperateTrajectory"),
   ...
) {
   model = match.arg(model)
   orig_args = x$`_function.params`
   orig_args$model = model

   more_args <- list(...)
   for(arg in names(more_args)) {
      orig_args[[arg]] <- more_args[[arg]]
   }
   #c(mean, more.args[!names(more.args) %in% names(mean)])

   do.call(ena, orig_args)
}

#' Title
#'
#' @param x [TBD]
#' @param by [TBD]
#' @param ... [TBD]
#'
#' @return [TBD]
#' @export
project_in <- function(x, by = NULL, ...) {
   if(is.null(by)) {
      stop("A second parameter (ena.set or rotation.set) is required")
   }

   rotation.set <- NULL
   if(is(by, "ena.set")) {
      rotation.set <- by$rotation
   } else if(is(by, "ena.rotation.set")) {
      rotation.set <- by
   }

   if(!identical(x$rotation$adjacency.key, rotation.set$adjacency.key)) {
      stop("Rotation sets must have identical adjacency keys")
   }

   x$rotation.matrix <- rotation.set$rotation.matrix
   x$rotation$rotation.matrix <- rotation.set$rotation.matrix
   x$rotation$nodes <- rotation.set$nodes;
   x$rotation$eigenvalues <- rotation.set$eigenvalues

   points <- as.matrix(x$model$points.for.projection) %*% as.matrix(x$rotation.matrix)
   points.dt <- as.data.table(points)
   for (i in seq(ncol(points.dt))) {
    set(points.dt, j = i, value = as.ena.dimension(points.dt[[i]]))
   }
   if(grepl(x = x$model$model.type, pattern = "Trajectory")) {
    x$points <- cbind(x$trajectories, points.dt)
   } else {
    x$points <- cbind(x$meta.data, points.dt)
   }
   x$points <- as.ena.matrix(x$points, "ena.points")

   .return(x, invisible = T)
}

#' Title
#'
#' @param x [TBD]
#' @param on [TBD]
#'
#' @return [TBD]
#' @export
means_rotate <- function(x, on = NULL) {
   groupVar = NULL
   groups = NULL
   if(is.null(on)) {
      col_counts = as.numeric(x$model$raw.input[, lapply(.SD, function(s) {
                  length(unique(s))
               }),
               .SDcols = c(x$`_function.params`$units)
            ])
      groupVar = x$`_function.params`$units[order(col_counts) == 1]
      groups = levels(unique(x$model$raw.input[[groupVar]]))[1:2]
      # on_grps = list()
      # on_grps[[on]] = sapply(on_vals, function(v) {
      #    x$meta.data[[on]] == v
      # }, simplify = F)
   } else if(!is.null(names(on))) {
      groupVar = names(on)
      groups = on[[groupVar]]
   }

   if(is.null(groupVar) || is.null(groups)) {
      stop("Unable to determine groups for rotation.")
   }

   orig_args <- x$`_function.params`
   orig_args$groupVar = groupVar
   orig_args$groups = groups
   new_set <- do.call(ena, orig_args)
   new_set$plots <- x$plots
   invisible(new_set)
}

.return <- function(x, invisible = T, from_plot = F) {
   # browser()
   x$`_plot_op` = from_plot
# if() {
#       print(x$plots)
#    }

   if(invisible == T) {
      invisible(x)
   } else {
      return(x)
   }
}
"""