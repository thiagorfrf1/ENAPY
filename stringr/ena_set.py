ena_set="""
ena.set <- function(x) {
  newset = list()
  class(newset) <- c("ena.set", class(newset))
  x.is.set <- T
  if("ENAdata" %in% class(x)) {
    x <- list(enadata = x);
    x.is.set <- F
  }
  code.columns <- apply(x$enadata$adjacency.matrix, 2, paste, collapse = " & ")

  newset$connection.counts <- x$enadata$adjacency.vectors;
  colnames(newset$connection.counts) <- code.columns
  for (i in seq(ncol(newset$connection.counts))) {
    set(newset$connection.counts, j = i,
        value = as.ena.co.occurrence(newset$connection.counts[[i]]))
  }

  if (grepl(x = x$enadata$model, pattern = "Traj", ignore.case = T)) {
    newset$meta.data <- data.table::copy(x$enadata$trajectories$units)
    newset$meta.data[, ENA_UNIT := apply(x$enadata$trajectories$units, 1,
                        paste, collapse = ".")]

    newset$trajectories <- cbind(newset$meta.data, x$enadata$trajectories$step)
    for (i in seq(ncol(newset$trajectories))) {
      set(newset$trajectories, j = i,
          value = as.ena.metadata(newset$trajectories[[i]]))
    }
  }
  else {
    newset$meta.data <- x$enadata$metadata
  }

  if (!is.null(newset$meta.data) && ncol(newset$meta.data) > 0) {
    for (i in seq(ncol(newset$meta.data))) {
      set(newset$meta.data, j = i,
          value = as.ena.metadata(newset$meta.data[[i]]))
    }
  }

  if (x.is.set) {
    newset$line.weights <- as.data.table(cbind(
                              x$enadata$metadata, x$line.weights))

    newset$points <- cbind(x$enadata$metadata, x$points.rotated)
    newset$rotation.matrix <- x$rotation.set$rotation
  }

  newset$connection.counts <- cbind(newset$meta.data, newset$connection.counts)
  class(newset$connection.counts) <- c("ena.connections",
                                        class(newset$connection.counts))

  newset$model <- list(
    model.type = x$enadata$model,
    raw.input = x$enadata$raw,
    row.connection.counts = x$enadata$accumulated.adjacency.vectors[,
            unique(names(x$enadata$accumulated.adjacency.vectors)), with = F],
    unit.labels = x$enadata$unit.names
  )

  # if(quote(x$enadata$function.params$weight.by) != "binary") {
  #   newset$model$unweighted.connection.counts <- x$enadata$adjacency.vectors.raw
  #   class(newset$model$unweighted.connection.counts) <- c("ena.connections",
  #                             class(newset$model$unweighted.connection.counts))
  #   are.codes <- find_code_cols(newset$model$unweighted.connection.counts)
  #   for (i in seq(are.codes)) {
  #     if (are.codes[i]) {
  #       set(newset$model$unweighted.connection.counts, j = i,
  #         value = as.ena.co.occurrence(
  #           newset$model$unweighted.connection.counts[[i]]
  #         )
  #       )
  #     } else {
  #       set(newset$model$unweighted.connection.counts, j = i,
  #         value = as.ena.metadata(
  #           newset$model$unweighted.connection.counts[[i]]
  #         )
  #       )
  #     }
  #   }
  # }

  cols <- grep("adjacency.code", colnames(newset$model$row.connection.counts))
  colnames(newset$model$row.connection.counts)[cols] <- code.columns
  for(i in cols) {
    set(newset$model$row.connection.counts, j = i,
        value = as.ena.co.occurrence(newset$model$row.connection.counts[[i]]))
  }
  for (i in which(colnames(newset$model$row.connection.counts)
      %in% colnames(newset$meta.data))
  ) {
    set(newset$model$row.connection.counts, j = i,
          value = as.ena.metadata(newset$model$row.connection.counts[[i]]))
  }
  for (i in which(colnames(newset$model$row.connection.counts) %in%
        x$enadata$codes)
  ) {
    set(newset$model$row.connection.counts, j = i,
          value = as.ena.code(newset$model$row.connection.counts[[i]]))
  }
  class(newset$model$row.connection.counts) <- c("row.connections",
                                      class(newset$model$row.connection.counts))

  if (x.is.set) {
    newset$model$centroids <- x$centroids
    newset$model$correlations <- x$correlations
    newset$model$function.call <- x$function.call
    newset$model$function.params <- x$function.params
    newset$model$points.for.projection <- cbind(x$enadata$metadata,
                                               x$points.normed.centered)
    newset$model$variance <- x$variance
    names(newset$model$variance) <- colnames(newset$rotation.matrix)
  }

  newset$rotation <- list(
    adjacency.key = as.data.table(x$enadata$adjacency.matrix),
    codes = x$enadata$codes
  )
  class(newset$rotation) <- c("ena.rotation.set", class(newset$rotation))

  for (i in seq(ncol(newset$rotation$adjacency.key))) {
    set(newset$rotation$adjacency.key, j = i,
          value = as.ena.codes(newset$rotation$adjacency.key[[i]]))
  }

  if(x.is.set) {
    newset$rotation$eigenvalues = x$rotation.set$eigenvalues
    newset$rotation$nodes = x$node.positions
    newset$rotation$rotation.matrix = x$rotation.set$rotation
  }

  newset$`_function.call` <- sys.calls()[[1]]
  back.frame <- sapply(sys.frames(), function(f) {
                                  "window.size.back" %in% ls(envir = f) })
  if (any(back.frame)) {
    call.frame <- sys.frame(which(back.frame))
    newset$`_function.params` <- mget(ls(envir = call.frame),
                                                  envir = call.frame)
  } else {
    newset$`_function.params` <- list()
  }

  return(newset);
}
"""