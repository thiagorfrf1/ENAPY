ena_plotter = """
#####
#' @title Wrapper to generate plots of units, groups, and networks
#'
#' @description Plots individual units, all units, groups of units, networks, and network subtractions
#'
#' @details This function includes options to plots individual units, all units,
#' groups of units, networks, and network subtractions, given an ena.set objects. Plots are stored
#' on the supplied ena.set object.
#'
#'
#' @param set an ena.set object
#' @param groupVar vector, character, of column name containing group identifiers.
#' @param groups vector, character, of values of groupVar column you wish to plot. Maxium of two groups allowed.
#' @param points logical, TRUE will plot points (default: FALSE)
#' @param mean logical, TRUE will plot the mean position of the groups defined in the groups argument (default: FALSE)
#' @param network logical, TRUE will plot networks (default: TRUE)
#' @param networkMultiplier numeric, scaling factor for non-subtracted networks (default: 1)
#' @param subtractionMultiplier numeric, scaling factor for subtracted networks (default: 1)
#' @param unit vector, character, name of a single unit to plot
#' @param print.plots logical, TRUE will show plots in the Viewer (default: FALSE)
#' @param ... Additional parameters passed to set creation and plotting functions
#' @export
#' @return ena.set object
#####
ena.plotter = function(
  set,
  groupVar = NULL,
  groups = NULL,
  points = FALSE,
  mean = FALSE,
  network = TRUE,
  networkMultiplier = 1,
  subtractionMultiplier = 1,
  unit = NULL,
  print.plots = F,
  ...
) {
  data = set$connection.counts;

  # set$plots[[length(set$plots)]] <- plot
  # plot <- set$plots[[length(set$plots)]]
  if(is.null(unit) == FALSE) {
    plot = ena.plot(enaset = set,title = unit)

    if(any(set$points$ENA_UNIT == unit) == FALSE){
      stop("Unit does not exist!")
    }

    point.row = set$points$ENA_UNIT == unit
    point = as.matrix(set$points)[point.row,]
    point.lw = as.matrix(set$line.weights)[point.row,]*networkMultiplier

    plot = ena.plot.points(enaplot = plot,points = point, colors = "black")
    plot = ena.plot.network(enaplot = plot, network = point.lw, colors = "black")

    set$plots[[length(set$plots) + 1]] <- plot

    if(print.plots == TRUE) {
      print(set$plots)
    }

    return(set)
  }

  if(is.null(groupVar) == TRUE) {
    plot = ena.plot(enaset = set, title = "All Units")

    if(network == TRUE) {
      lineweights = as.matrix(set$line.weights)
      mean.lineweights = colMeans(lineweights) * networkMultiplier

      plot = ena.plot.network(plot, network = mean.lineweights, colors = "black")
    }

    if(points == TRUE) {
      points.for.plot = as.matrix(set$points)

      plot = ena.plot.points(enaplot = plot,points = points.for.plot,colors = "black")
    }

    if(mean == TRUE) {
      points.for.plot = as.matrix(set$points)

      plot = ena.plot.group(plot, points.for.plot, colors = "black", labels = "Mean",confidence.interval = "box")
    }

    else if(TRUE %in% c(network,points, mean) == FALSE) {
      stop("You must set at least one of points, mean, or network to TRUE to obtain a plot.")
    }

    set$plots[[length(set$plots) + 1]] <- plot

    if(print.plots == TRUE) {
      print(set$plots)
    }

    return(set)
  }
  else if(is.null(groups) == TRUE) {
    unique.groups = unique(data[[groupVar]])

    if(length(unique.groups) == 1){
      warning("No groups specified and group variable only contains one unique value. Generating plot for one group.")

      group = unique.groups

      group.rows = set$points[[groupVar]] == group
      g.plot = ena.plot(enaset = set, title = group)

      if(network == TRUE) {
        g.lw = as.matrix(set$line.weights)[group.rows, , drop = FALSE]
        g.mean.lw = colMeans(g.lw) * networkMultiplier
        g.plot = ena.plot.network(g.plot, network = g.mean.lw, colors = "black")
      }

      if(points == TRUE) {
        g.points.for.plot = as.matrix(set$points)[group.rows, , drop = FALSE]
        g.plot = ena.plot.points(enaplot = g.plot,points = g.points.for.plot,colors = "black")
      }

      if(mean == TRUE) {
        g.points.for.plot = as.matrix(set$points)[group.rows, , drop = FALSE]
        g.plot = ena.plot.group(g.plot, g.points.for.plot, colors = "black", labels = group,confidence.interval = "box")
      }

      else if(TRUE %in% c(network,points, mean) == FALSE) {
        stop("You must set at least one of points, mean, or network to TRUE to obtain a plot.")
      }
      set$plots[[length(set$plots) + 1]] <- g.plot

      if(print.plots == TRUE) {
        print(set$plots)
      }

      return(set)
    }
    else {
      group1 = unique.groups[1]
      group2 = unique.groups[2]

      warning(paste0("No groups specified. Generating plots of first two unique values of group variable: ",group1," and ",group2))

      set = ena.plot.subtraction(set = set,
               groupVar = groupVar,
               group1 = group1,
               group2 = group2,
               points = points,
               mean = mean,
               network = network,
               networkMultiplier = networkMultiplier,
               subtractionMultiplier = subtractionMultiplier)


      if(print.plots == TRUE) {
        print(set$plots)
      }

       return(set)
    }
  }
  else if(length(groups) == 1) {
    group = groups

    if(any(data[[groupVar]] == group) == FALSE){
      stop("Group column does not contain group1 value!")
    }

    group.rows = set$points[[groupVar]] == group
    g.plot = ena.plot(enaset = set, title = group)

    if(network == TRUE) {
      g.lw = as.matrix(set$line.weights)[group.rows, , drop = FALSE]
      g.mean.lw = colMeans(g.lw) * networkMultiplier

      g.plot = ena.plot.network(g.plot, network = g.mean.lw, colors = "black")
    }

    if(points == TRUE) {
      g.points.for.plot = as.matrix(set$points)[group.rows, , drop = FALSE]
      g.plot = ena.plot.points(enaplot = g.plot,points = g.points.for.plot,colors = "black")
    }

    if(mean == TRUE) {
      g.points.for.plot = as.matrix(set$points)[group.rows, , drop = FALSE]
      g.plot = ena.plot.group(g.plot, g.points.for.plot, colors = "black", labels = group,confidence.interval = "box")
    }

    else if(TRUE %in% c(network,points, mean) == FALSE) {
      stop("You must set at least one of points, mean, or network to TRUE to obtain a plot.")
    }
    set$plots[[length(set$plots) + 1]] <- g.plot

    if(print.plots == TRUE) {
      print(set$plots)
    }

    return(set)
  }
  else if (length(groups) >= 2) {
    if (length(groups) > 2) {
      warning(paste0("More than two groups specified. Plotting the first two groups: ", groups))
    }

    groups.missing = groups[which(!groups %in% data[[groupVar]])]
    if(length(groups.missing) > 0) {
      stop(paste0("Group column does not contain group value(s): ", groups[groups.missing]))
    }

    set = ena.plot.subtraction(
      set = set,
      groupVar = groupVar,
      group1 = groups[1],
      group2 = groups[2],
      points = points,
      mean = mean,
      network = network,
      networkMultiplier = networkMultiplier,
      subtractionMultiplier = subtractionMultiplier,
      ...
    )

    if(print.plots == TRUE) {
      print(set$plots)
    }

    return(set)
  }
}

"""