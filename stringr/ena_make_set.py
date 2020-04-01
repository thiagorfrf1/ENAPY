ena_make_set = """
##
#' @title Generate ENA Set
#'
#' @description Generates an ENA model by constructing a dimensional reduction of adjacency (co-occurrence) vectors in an ENA data object
#'
#' @details This function generates an ENAset object from an ENAdata object. Takes
#' the adjacency (co-occurrence) vectors from enadata, computes a dimensional
#' reduction (projection), and calculates node positions in the projected ENA
#' space. Returns location of the units in the projected space, as well as
#' locations for node positions, and normalized adjacency (co-occurrence) vectors
#' to construct network graphs
#'
#' @export
#'
#' @param enadata \code{\link{ENAdata}} that will be used to generate an ENA model
#' @param dimensions The number of dimensions to include in the dimensional reduction
#' @param norm.by A function to be used to normalize adjacency (co-occurrence) vectors before computing the dimensional reduction, default: sphere_norm_c()
#' @param rotation.by	A function to be used to compute the dimensional reduction, default: ena.svd()
#' @param rotation.params (optional) A character vector containing additional parameters for the function in rotation.by, if needed
#' @param rotation.set A previously-constructed  ENARotationSet object to use for the dimensional reduction
#' @param endpoints.only A logical variable which determines whether to only show endpoints for trajectory models
#' @param node.position.method A function to be used to determine node positions based on the dimensional reduction, default: lws.position.es()
#' @param as.list R6 objects will be deprecated, but if this is TRUE, the original R6 object will be returned, otherwise a list with class `ena.set`
#' @param ... additional parameters addressed in inner function
#'
#' @examples
#' data(RS.data)
#'
#' codeNames = c('Data','Technical.Constraints','Performance.Parameters',
#'   'Client.and.Consultant.Requests','Design.Reasoning','Collaboration');
#'
#' accum = ena.accumulate.data(
#'   units = RS.data[,c("UserName","Condition")],
#'   conversation = RS.data[,c("Condition","GroupName")],
#'   metadata = RS.data[,c("CONFIDENCE.Change","CONFIDENCE.Pre","CONFIDENCE.Post")],
#'   codes = RS.data[,codeNames],
#'   window.size.back = 4
#' )
#'
#' set = ena.make.set(
#'   enadata = accum
#' )
#'
#' set.means.rotated = ena.make.set(
#'   enadata = accum,
#'   rotation.by = ena.rotate.by.mean,
#'   rotation.params = list(
#'       accum$meta.data$Condition=="FirstGame",
#'       accum$meta.data$Condition=="SecondGame"
#'   )
#' )
#'
#' @seealso \code{\link{ena.accumulate.data}}, \code{\link{ENAset}}
#'
#' @return \code{\link{ENAset}} class object that can be further processed for analysis or plotting
##
ena.make.set <- function(
  enadata,
  dimensions = 2,
  norm.by = fun_sphere_norm,
  rotation.by = ena.svd,
  rotation.params = NULL,
  rotation.set = NULL,
  endpoints.only = T,
  node.position.method = lws.positions.sq,
  as.list = TRUE,
  ...
) {
  if (as.list == F) {
    warning(paste0("Usage of ENAdata and ENAset objects will be deprecated ",
      "and potentially removed altogether in future versions."))

    if (!is(enadata, "ENAdata")) {
      stop(paste0("Use of ena.make.set with as.list=FALSE requires `enadata` ",
        "be an ENAdata object. Re-run the accumulation with as.list=FALSE"))
    }

    set <- ENAset$new(
      enadata = enadata,
      dimensions = dimensions,
      rotation.by = ifelse(
        identical(rotation.by, ena.svd),
        ena.svd.R6,
        rotation.by
      ),
      rotation.params = rotation.params,
      rotation.set = rotation.set,
      norm.by = norm.by,
      node.position.method = ifelse(
        identical(node.position.method, lws.positions.sq),
        lws.positions.sq.R6,
        node.position.method
      ),
      endpoints.only = endpoints.only,
      ...
    )
    return(set$process());
  }
  else {
    if ("ENAdata" %in% class(enadata)) {
      warning(paste0("Usage of ENAdata objects will be deprecated and ",
        "potentially removed altogether in future versions. See ",
        "ena.accumulate.data() or ena.set()."))

      enadata <- ena.set(enadata)
    }

    ###
    # Convert the string vector of code names to their corresponding
    # co-occurence names
    #####
      code_columns <- svector_to_ut(enadata$rotation$codes)

    ###
    # Normalize the raw data using self$function.params$norm.by,
    # which defaults to calling rENA::dont_sphere_norm_c
    #####
      line.weights <- norm.by(as.matrix(enadata$connection.counts))
      colnames(line.weights) <- code_columns

      line.weights.dt <- as.data.table(line.weights)
      for (i in seq(ncol(line.weights.dt)))
        set(line.weights.dt, j = i,
            value = as.ena.co.occurrence(line.weights.dt[[i]]))

      enadata$line.weights <- cbind(enadata$meta.data, line.weights.dt)
      class(enadata$line.weights) <- c("ena.line.weights",
                                      class(enadata$line.weights))
    #####

    ###
    # Center the normed data
    #####
      points.for.projection <- center_data_c(line.weights)
      colnames(points.for.projection) <- code_columns;
      enadata$model$points.for.projection = as.data.table(points.for.projection)
      for (i in seq(ncol(enadata$model$points.for.projection))) {
        set(
          enadata$model$points.for.projection,
          j = i,
          value = as.ena.co.occurrence(enadata$model$points.for.projection[[i]])
        )
      }
      enadata$model$points.for.projection <- as.ena.matrix(cbind(
        enadata$meta.data,
        enadata$model$points.for.projection
      ), "ena.points")
    #####

    ###

    ###
    # Generate and Assign the rotation set
    #####
      if (!is.null(rotation.by) && is.null(rotation.set)) {
        rotation <- do.call(rotation.by, list(enadata, rotation.params))

        enadata$rotation.matrix <- as.data.table(rotation$rotation, keep.rownames = "codes")
        for (i in seq(ncol(enadata$rotation.matrix))) {
          if(i == 1) {
            set(enadata$rotation.matrix,
                j = i, value = as.ena.metadata(enadata$rotation.matrix[[i]])
            )
          } else {
            set(enadata$rotation.matrix,
                j = i, value = as.ena.dimension(enadata$rotation.matrix[[i]])
            )
          }
        }
        class(enadata$rotation.matrix) <- c("ena.rotation.matrix", class(enadata$rotation.matrix))

        enadata$rotation$rotation.matrix <- enadata$rotation.matrix
        enadata$rotation$eigenvalues <- rotation$eigenvalues;
      }
      else if (!is.null(rotation.set)) {
        if (is(rotation.set, "ena.rotation.set")) {
          enadata$rotation.matrix <- rotation.set$rotation.matrix
          enadata$rotation$rotation.matrix <- rotation.set$rotation.matrix
          enadata$rotation$nodes <- rotation.set$nodes;
          enadata$rotation$eigenvalues <- rotation.set$eigenvalues
        } else {
          stop("Supplied rotation.set is not an instance of ENARotationSet")
        }
      }
      else {
        stop("Unable to find or create a rotation set")
      }
    #####

    ###
    # Generate the rotated points
    #####
      if (!is.null(enadata$rotation.matrix)) {
        points <- points.for.projection %*% as.matrix(enadata$rotation.matrix)
        points.dt <- as.data.table(points)
        for (i in seq(ncol(points.dt))) {
          set(points.dt, j = i, value = as.ena.dimension(points.dt[[i]]))
        }
        if(grepl(x = enadata$model$model.type, pattern = "Trajectory")) {
          enadata$points <- cbind(enadata$trajectories, points.dt)
        }
        else {
          enadata$points <- cbind(enadata$meta.data, points.dt)
        }
        enadata$points <- as.ena.matrix(enadata$points, "ena.points")
      }
      else {
        stop(paste0("There is no rotation matrix, if you supplied a custom ",
          "rotation.set, be sure it contains a rotation.matrix"))
      }
    #####

    ###
    # Calculate node positions
    #  - The supplied methoed is responsible is expected to return a list
    #    with two keys, "node.positions" and "centroids"
    #####
      if (exists("rotation") && !is.null(rotation) && is.null(rotation.set)) {
        positions <- node.position.method(enadata)

        if (all(names(positions) %in% c("node.positions", "centroids"))) {
          enadata$rotation$nodes <- as.data.table(positions$node.positions)
          colnames(enadata$rotation$nodes) <- colnames(points)
          rownames(enadata$rotation$nodes) <- enadata$rotation$codes

          for (i in seq(ncol(enadata$rotation$nodes))) {
            set(enadata$rotation$nodes, j = i,
                  value = as.ena.dimension(enadata$rotation$nodes[[i]]))
          }
          enadata$rotation$nodes <- data.table(
            code = structure(enadata$rotation$codes,
                      class = c("code", class(enadata$rotation$codes))),
            enadata$rotation$nodes
          )
          class(enadata$rotation$nodes) = c("ena.nodes",
                                            class(enadata$rotation$nodes))

          enadata$model$centroids <- as.data.table(positions$centroids)
          for (i in seq(ncol(enadata$model$centroids))) {
            set(enadata$model$centroids, j = i,
              value = as.ena.dimension(enadata$model$centroids[[i]])
            )
          }
          colnames(enadata$model$centroids) <- colnames(as.matrix(enadata$rotation.matrix))
          enadata$model$centroids = cbind(
            data.table(unit = enadata$model$unit.labels),
            enadata$model$centroids
          )
          set(enadata$model$centroids, j = 1L,
            value = as.ena.metadata(enadata$model$centroids[[1L]])
          )
          enadata$model$centroids <- as.ena.matrix(enadata$model$centroids)
        }
        else {
          stop(paste0("The node position method didn't return back the ",
            "expected objects:\n",
            "\tExpected: c('node.positions','centroids')\n",
            "\tReceived: ", names(positions), sep = ""))
        }
      } else if (!is.null(rotation.set)) {
        enadata$rotation$nodes <- rotation.set$nodes
      }

      if (is.null(enadata$rotation$nodes)) {
        stop("Unable to determine the node positions either by calculating
                    them using `node.position.method` or using a supplied
                    `rotation.set`")
      }
    #####

    ###
    # Variance
    #####
      var_rot_data <- var(points)
      diagonal_variance <- as.vector(diag(var_rot_data))
      enadata$model$variance <- diagonal_variance / sum(diagonal_variance)
      names(enadata$model$variance) <- colnames(enadata$rotation$rotation.matrix)[-1]
    #####

    enadata$plots <- list() #default = ena.plot(enadata, ...))
    # class(enadata$model$plot) <- c("ena.plot", class(enadata$model$plot))

    enadata$`_function.params`$norm.by <- norm.by

    return(enadata)
  }
}
"""