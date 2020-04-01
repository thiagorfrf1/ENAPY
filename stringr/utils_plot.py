utils_plot = """

#####
#' Plot an ena.set object
#'
#' @param x ena.set to plot
#' @param y ignored.
#' @param ... Additional parameters passed along to ena.plot functions
#'
#' @examples
#' library(magrittr)
#'
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
#' plot(set) %>%
#'   add_points(Condition$FirstGame, colors = "blue", with.mean = TRUE) %>%
#'   add_points(Condition$SecondGame, colors = "red", with.mean = TRUE)
#'
#' plot(set) %>%
#'   add_network(Condition$FirstGame - Condition$SecondGame)
#'
#' @return ena.plot.object
#' @export
#####
plot.ena.set <- function(x, y, ...) {
  p = ena.plot(x, ...)
  # p
  # p$enaset = NULL
  x$plots[[length(x$plots) + 1]] = p
  args = list(...)
  if(!is.null(args$title)) {
    names(x$plots)[length(x$plots)] = args$title
  }

  .return(x, from_plot = T, invisible = F)
}

#' Plot points on an ena.plot
#'
#' @param x ena.plot to add point on
#' @param wh which points to plot
#' @param ... additional parameters to pass along
#' @param name name to give the plot
#' @param mean include a mean point for the provided points
#' @param colors colors for plotted points
#'
#' @return ena.plot.object
#' @export
add_points <- function(
  x,
  wh = NULL, ...,
  name = "plot",
  mean = NULL,
  colors = NULL
) {
  set <- x
  plot <- set$plots[[length(set$plots)]]
  more.args <- list(...)

  wh_subbed <- as.character(substitute(wh))
  if (!is.null(wh_subbed) && length(wh_subbed) > 0) {
    if (length(wh_subbed) > 1 && wh_subbed[[2]] %in% colnames(set$points)) {
      cc <- call(wh_subbed[[1]], set$points, wh_subbed[[2]])
      part1 <- eval(cc)

      name <- paste(wh_subbed[-1], collapse = "$")
      if(grepl(set$model$model.type, pattern="Trajectory")) {
        points <- set$points[part1 == wh_subbed[[3]], ]
        more.args$points = points[, .SD[nrow(.SD)], by = ENA_UNIT]
      }
      else {
        more.args$points = points <- set$points[part1 == wh_subbed[[3]], ]
      }

      if(is.null(colors)) {
        colors = plot$palette[length(plot$plotted$points) + 1]
      }
    }
    else if (length(wh_subbed) == 1 && wh_subbed[[1]] %in% colnames(set$points)) {
      more.args$points = points = set$points
      if(is.null(colors)) {
        colors <- plot$palette[as.numeric(as.factor(set$points[[wh_subbed]])) + length(plot$plotted$points)]
      }
      else {
        colors <- colors[as.numeric(as.factor(set$points[[wh_subbed]]))]
      }
    }
    else {
      more.args$points = points <- wh
      colors = ifelse(is.null(colors), plot$palette[length(plot$plotted$points) + 1], colors)
    }
  }
  else {
    more.args$points = points = set$points
    name <- "all.points"
    colors = ifelse(is.null(colors), plot$palette[length(plot$plotted$points) + 1], colors)
  }

  more.args$enaplot = plot
  more.args$legend.name = name
  if(!is.null(colors)) {
    more.args$colors = colors
  } else {
    more.args$colors = plot$palette[length(plot$plotted$points) + 1]
  }
  plot <- do.call(ena.plot.points, more.args)

  for(color in unique(more.args$colors)) {
    plot$plotted$points[[length(plot$plotted$points) + 1]] <- list(
      data = more.args$points[color == more.args$colors,],
      color = color
    )
    if(!is.null(name)) {
      names(plot$plotted$points)[length(plot$plotted$points)] = name
    }
  }

  if(!is.null(mean) && (is.list(mean) || mean == T)) {
    if (is.list(mean)) {
      more.args <- c(mean, more.args[!names(more.args) %in% names(mean)])
    }
    more.args$enaplot <- plot
    more.args$points <- points
    more.args$labels <- name

    plot <- do.call(ena.plot.group, more.args)
  }

  set$plots[[length(set$plots)]] <- plot
  invisible(set)
}

#' Plot a trajectory on an ena.plot
#'
#' @param x ena.plot object to plot on
#' @param wh which points to plot as the trajectory
#' @param name Name, as a character vector, to give the plot
#' @param ... additional parameters to pass along
#'
#' @return ena.plot.object
#' @export
add_trajectory <- function(x, wh = NULL, ..., name = "plot") {
  set <- x
  # plot <- set$model$plot
  plot <- set$plots[[length(set$plots)]]

  subbed <- substitute(wh)
  args_list <- as.character(subbed)
  points <- set$points

  if (!is.null(args_list) && !is.null(subbed)) {
    if (length(args_list) > 1) {
      wh_subbed <- as.character(substitute(wh))
      cc <- call(wh_subbed[[1]], set$points, wh_subbed[[2]])
      part1 <- eval(cc)
      points <- set$points[part1 == wh_subbed[[3]], ]
      by <- "ENA_UNIT"
    }
    else {
      by <- args_list[[1]]
    }
  }
  else {
    by <- "ENA_UNIT"
  }
  plot <- ena.plot.trajectory(plot, points = points, by = by)

  # set$model$plot <- plot
  set$plots[[length(x$plots)]] <- plot
  invisible(set)
}

#' Add a group mean to an ena.plot
#'
#' @param x ena.plot object to plot on
#' @param wh which points to plot as the trajectory
#' @param ... additional parameters to pass along
#'
#' @return ena.plot.object
#' @export
add_group <- function(x, wh = NULL, ...) {
  set <- x
  # plot <- set$model$plot
  plot <- set$plots[[length(set$plots)]]

  arg_list <- list(...)
  wh.clean <- substitute(wh)

  if (
    identical(as.character(wh.clean), "wh.clean") ||
    identical(as.character(wh.clean), "y")
  ) {
    wh.clean <- wh;
  }

  more_args = list(...)
  more_args$enaplot <- plot
  if(is.null(more_args$color)) {
    more_args$colors = plot$palette[length(plot$plotted$points) + 1]
  }

  if (is.null(wh.clean)) {
    plot <- do.call(ena.plot.group, more_args)
  }
  else {
    parts <- as.character(wh.clean)

    if (parts[2] %in% colnames(set$line.weights)) {
      label <- parts[3]
      group.rows <- set$points[set$points[[parts[2]]] == parts[3], ]
      if(nrow(group.rows) > 0) {
        group.means <- colMeans(group.rows)

        more_args$points <- group.means
        more_args$labels <- label
        plot <- do.call(ena.plot.group, more_args)
      }
      else {
        warning("No points in the group")
      }
    }
    else {
      warning("Unable to plot group")
    }
  }

  plot$plotted$means[[length(plot$plotted$means) + 1]] = list(
    data = more_args$points,
    color = more_args$colors
  )

  set$plots[[length(set$plots)]] <- plot
  invisible(set)
}

#' Add a network to an ENA plot
#'
#' @param x ena.plot object to plot wtih
#' @param wh network to plot
#' @param with.mean Logical value, if TRUE plots the mean for the points in the network
#' @param ... Additional parametesr to pass along
#'
#' @return ena.plot.object
#' @export
add_network <- function(x, wh = NULL, ..., with.mean = F) {
  set <- x
  # plot <- set$model$plot
  plot <- set$plots[[length(set$plots)]]

  wh.clean <- substitute(wh)
  arg_list <- list(...)

  if(is.null(wh.clean)) { #, "ena.points")) {
    plot <- ena.plot.network(
      plot,
      network = colMeans(set$line.weights),
      points = set$rotation$nodes[, 1:2],
      ...
    )

    if (with.mean) {
      set <- add_group(set, points = set$points, ...)
      plot <- set$plots[[length(set$plots)]]
    }
  }
  else {
    parts <- as.character(wh.clean)

    if (length(wh.clean) > 1 && is.call(wh.clean[[2]])) {
      means <- sapply(c(wh.clean[[2]], wh.clean[[3]]), function(y) {
        parts <- as.character(y)

        if(with.mean) {
          set <- add_group(set, y,
                colors = plot$palette[length(attr(plot, "means")) + 1], ...)
          plot <- set$plots[[length(set$plots)]]
        }

        colMeans(set$line.weights[set$line.weights[[parts[2]]] == parts[3], ])
      })

      group.means <- means[, 1] - means[, 2]
    }
    else {
      if (parts[2] %in% colnames(set$line.weights)) {
        group.means <- colMeans(
          as.matrix(set$line.weights[set$line.weights[[parts[2]]] == parts[3], ])
        )

        if (with.mean) {
          set <- add_group(set, wh.clean, ...)
          plot <- set$plots[[length(set$plots)]]
        }
      }
      else {
        wgts <- get(as.character(wh.clean), envir = parent.frame())
        group.means <- colMeans(wgts)
        if (with.mean) warning("Not able to determine mean automatically")
      }
    }

    plot <- ena.plot.network(plot,
          network = group.means,
          node.positions = as.matrix(set$rotation$nodes)[, 1:2], ...)
  }

  # set$model$plot <- plot
  set$plots[[length(set$plots)]] <- plot
  invisible(set)
}

#' Title
#'
#' @param x [TBD]
#' @param ... [TBD]
#'
#' @return [TBD]
#' @export
add_nodes <- function(x, ...) {
  set <- x
  plot <- set$plots[[length(set$plots)]]

  nodes <- set$rotation$nodes
  plot <- ena.plot.points(plot,
            points = as.matrix(nodes),
            texts = as.character(nodes$code),
            ...
          )

  plot$plotted$networks[[length(plot$plotted$networks) + 1]] <- list(
    nodes = nodes,
    data = NULL,
    color = NULL
  )
  set$plots[[length(set$plots)]] <- plot
  invisible(set)
}

#' Title
#'
#' @param x [TBD]
#'
#' @return [TBD]
#' @export
with_means <- function(x) {
  set <- x
  # plot <- set$model$plot
  plot <- set$plots[[length(set$plots)]]

  for(point_group in plot$plotted$points) {
    plot <- ena.plot.group(plot, point_group$data, colors = point_group$color[1])

    plot$plotted$means[[length(plot$plotted$means) + 1]] <- list(
      data = colMeans(point_group$data),
      color = point_group$color[1]
    )
  }

  # set$model$plot <- plot
  set$plots[[length(set$plots)]] <- plot
  invisible(set)
}

#' Title
#'
#' @param x [TBD]
#' @param ... [TBD]
#' @param by [TBD]
#' @param add_jitter [TBD]
#' @param frame [TBD]
#' @param transition [TBD]
#' @param easing [TBD]
#'
#' @return [TBD]
#' @export
with_trajectory <- function(
  x, ...,
  by = x$`_function.params`$conversation[1],
  add_jitter = TRUE,
  frame = 1100,
  transition = 1000,
  easing = "circle-in-out"
) {
  set <- x
  if(!grepl(x = set$model$model.type, pattern = "Trajectory")) {
    stop(paste0("Unable to plot trajectories on model of type: ", set$model$model.type))
  }
  plot <- set$plots[[length(set$plots)]]

  args = list(...)


  all_steps_w_zero <- data.table(rbind(
    rep(0, length(by)),
    expand.grid(
      sapply(by, function(b) sort(unique(set$points[[b]]))),
      stringsAsFactors = F
    )
  ))
  colnames(all_steps_w_zero) <- by
  point_group_names <- seq(plot$plotted$points)
  points_cleaned <- lapply(point_group_names, function(n) {
    prepare_trajectory_data(
      points = plot$plotted$points[[n]]$data,
      by = by,
      units = plot$plotted$points[[n]]$data,
      units_by = set$`_function.params`$units,
      steps = all_steps_w_zero
    )
  })
  names(points_cleaned) <- sapply(plot$plotted$points, "[[", "color")
  points_cleaned <- rbindlist(points_cleaned, idcol = "color")

  meta_data = unique(set$meta.data)
  setkey(points_cleaned, ENA_UNIT)
  setkey(meta_data, ENA_UNIT)
  points_cleaned = meta_data[points_cleaned]
  setkeyv(points_cleaned, by)

  size = ifelse(is.null(args$size), 10, args$size)
  opacity = ifelse(is.null(args$opacity), 1, args$opacity)

  dims = as.matrix(points_cleaned[, find_dimension_cols(points_cleaned), with = F])[, 1:2]
  if(add_jitter) {
    dims[, 1] = jitter(dims[, 1])
    dims[, 2] = jitter(dims[, 2])
  }

  if(is.null(args$scale)) {
    max_abs = max(abs(dims))
    scale = c(-1*max_abs, max_abs)
  }
  else {
    scale = args$scale
  }

  ax <- list(
    range = scale, title = "",
    zeroline = TRUE, showline = FALSE,
    showticklabels = FALSE, showgrid = FALSE
  )

  #####
  ### Add to the plot
  #####
    thisPlot <- plotly::plot_ly(
        data = points_cleaned,
        x = dims[,1], y = dims[,2],
        text = ~ENA_UNIT,
        frame = as.formula(paste0("~", by)),
        type = 'scatter',
        mode = 'markers',
        marker = list(
          size = size,
          opacity = opacity,
          hoverinfo = "text",
          color = as.numeric(as.factor(points_cleaned[["color"]]))
        )
      ) %>%
      plotly::layout(
        xaxis = ax,
        yaxis = ax,
        showlegend = T
      ) %>%
      plotly::animation_opts(
        frame = frame,
        transition = transition,
        easing = easing,
        redraw = T
      )
  #####

  # set$model$plot <- plot
  set$plots[[length(set$plots) + 1]] <- thisPlot
  invisible(set)
}


#' Title
#'
#' @param x [TBD]
#' @param by [TBD]
#' @param rotation_matrix [TBD]
#' @param points [TBD]
#' @param units [TBD]
#' @param units_by [TBD]
#' @param steps [TBD]
#'
#' @return [TBD]
#' @export
prepare_trajectory_data <- function(
  x = NULL,
  by = x$`_function.params`$conversation[1],
  rotation_matrix = x$rotation.matrix,
  points = NULL,
  units = points,
  units_by = x$`_function.params`$units,
  steps = NULL
) {
  if(is(x, "ena.set")) {
    if(is.null(points))
      points <- x$points
    if(is.null(units))
      units <- x$trajectories #points[, find_meta_cols(points), with = FALSE]
  }

  unique_unit_values <- unique(units[, c(units_by, "ENA_UNIT"), with = FALSE])

  if(!is.null(rotation_matrix)) {
    rotation_matrix = as.matrix(rotation_matrix)
    full_data <- cbind(units, as.matrix(points) %*% rotation_matrix)
  } else {
    full_data <- cbind(units, as.matrix(points))
  }
  full_data <- full_data[, unique(names(full_data)), with = FALSE]

  if(is.null(steps)) {
    all_steps_w_zero <- data.table(rbind(
      rep(0, length(by)),
      expand.grid(
        sapply(by, function(b) sort(unique(units[[b]]))),
        stringsAsFactors = F
      )
    ))
    colnames(all_steps_w_zero) <- by
  } else {
    all_steps_w_zero <- steps
  }
  all_step_data <- CJ(all_steps_w_zero[[by]], unique_unit_values$ENA_UNIT)
  colnames(all_step_data) <- c(by, "ENA_UNIT")

  dimension_col_names = colnames(points)[
                          which(sapply(points, function(col) {
                            is(col, "ena.dimension")
                          }))
                        ]
  all_step_data[, c(dimension_col_names) := 0]
  all_step_data[[by]] = as.ena.metadata(all_step_data[[by]])
  all_step_data = merge(unique_unit_values, all_step_data, by = "ENA_UNIT")
  setkey(all_step_data, "ENA_UNIT")

  filled_data = all_step_data[ , {
      by_names = names(.BY)
      user_rows = sapply(1:length(by_names), function(n) {
          full_data[[by_names[n]]] == .BY[n]
      })
      existing_row = which(rowSums(user_rows * 1) == 2)
      if(length(existing_row) > 0) {
        full_data[existing_row, c(dimension_col_names), with = FALSE]
      } else {
        prev_row = tail(full_data[ENA_UNIT == .BY$ENA_UNIT & full_data[[by]] < .BY[[by]],], 1)
        if(nrow(prev_row) == 0) {
          data.table(matrix(rep(0, length(dimension_col_names)), nrow = 1, dimnames = list(NULL, c(dimension_col_names))))
        } else {
          prev_row[, c(dimension_col_names), with = FALSE]
        }
      }

  },  by = c("ENA_UNIT", by)]
  for(col in dimension_col_names) {
    set(filled_data, j = col, value = as.ena.dimension(filled_data[[col]]))
  }
  return(filled_data)
}


#' Title
#'
#' @param x [TBD]
#' @param wh [TBD]
#'
#' @return [TBD]
#' @export
clear <- function(x, wh = seq(x$plots)) {
  if(length(wh) > 0) {
    x$plots[[wh]] <- NULL
  }
  invisible(x)
}

#' Title
#'
#' @param x [TBD]
#' @param center Ignored.
#' @param scale [TBD]
#'
#' @return [TBD]
#' @export
scale.ena.set <- function(x, center = TRUE, scale = TRUE) {
  set <- x
  plot <- set$plots[[length(set$plots)]]
  browser()
  dims <- 1:2
  point_range <- range(sapply(plot$plotted$points, function(d) range(as.matrix(d$data)[,dims])))
  network_range <-range(sapply(plot$plotted$networks, function(n) range(as.matrix(n$nodes)[,dims])))

  scale_factor <- min(abs(network_range) / abs(point_range))

  for( points in plot$plotted$points) {
    dim_cols = colnames(points$data)[find_dimension_cols(points$data)]
    points$data[, c(dim_cols) := lapply(.SD, function(x) x * scale_factor), .SDcols = c(dim_cols)]
    more_args = list()
    more_args$enaplot <- plot
    more_args$points <- points$data
    more_args$colors <- points$color
    plot <- do.call(ena.plot.points, more_args)
  }
  for(means in plot$plotted$means) {
    more_args <- list()
    more_args$enaplot <- plot
    more_args$points <- means$data * scale_factor
    more_args$colors <- means$color
    plot <- do.call(ena.plot.group, more_args)
  }

  set$plots[[length(set$plots)]] <- plot

  invisible(set)
}

check_range <- function(x) {
  numbers <- as.numeric(sapply(x$plotted$points, function(p) max(as.matrix(p$data))))
  network <- as.numeric(sapply(x$plotted$network, function(p) max(as.matrix(p$nodes))))
  means <- as.numeric(sapply(x$plotted$means, function(p) max(as.matrix(p$data))))

  if(
    length(numbers) == 0 &&
    length(means) == 0
  ) {
    return(x)
  }

  curr_max = max(c(numbers, network, means))
  if(curr_max*1.2 > max(x$axes$y$range)) {
    this.max = curr_max * 1.2
    x$axes$x$range = c(-this.max, this.max)
    x$axes$y$range = c(-this.max, this.max)
    x$plot = plotly::layout(
      x$plot,
      xaxis = x$axes$x,
      yaxis = x$axes$y
    );
  } else if (curr_max < max(x$axes$y$range*0.5)) {
    this.max = curr_max * 1.2
    x$axes$x$range = c(-this.max, this.max)
    x$axes$y$range = c(-this.max, this.max)
    x$plot = plotly::layout(
      x$plot,
      xaxis = x$axes$x,
      yaxis = x$axes$y
    );
  }

  x
}

#' Title
#'
#' @param x [TBD]
#' @param ... [TBD]
#'
#' @return [TBD]
#' @export
show <- function(x, ...) {
   x$plots <- lapply(x$plots, check_range)
   print(x, ..., plot = T, set = F)
   invisible(x)
}
"""