ena_plot_subtraction = """
### plot subtraction ###

ena.plot.subtraction = function(
  set,
  groupVar = NULL,
  group1 = NULL,
  group2 = NULL,
  points = FALSE,
  mean = FALSE,
  network = TRUE,
  networkMultiplier = 1,
  subtractionMultiplier = 1,
  ...
) {
  group1.rows = set$points[[groupVar]] == group1
  group2.rows = set$points[[groupVar]] == group2

  g1.plot = ena.plot(enaset = set, title = group1)
  g2.plot = ena.plot(enaset = set, title = group2)
  sub.plot = ena.plot(enaset = set, title = paste0("Network Subtraction -- ",group1," vs ",group2))

  if(network == TRUE) {
    g1.lw = as.matrix(set$line.weights)[group1.rows,,drop=FALSE]
    g1.mean.lw = colMeans(g1.lw) * networkMultiplier

    g2.lw = as.matrix(set$line.weights)[group2.rows,,drop=FALSE]
    g2.mean.lw = colMeans(g2.lw) * networkMultiplier

    sub = (g1.mean.lw - g2.mean.lw) * subtractionMultiplier

    g1.plot = ena.plot.network(g1.plot, network = g1.mean.lw, colors = "blue")
    g2.plot = ena.plot.network(g2.plot, network = g2.mean.lw, colors = "red")
    sub.plot = ena.plot.network(sub.plot, network = sub)
  }

  if(points == TRUE) {
    g1.points.for.plot = as.matrix(set$points)[group1.rows,,drop=FALSE]
    g2.points.for.plot = as.matrix(set$points)[group2.rows,,drop=FALSE]

    g1.plot = ena.plot.points(enaplot = g1.plot, points = g1.points.for.plot, colors = "blue")
    g2.plot = ena.plot.points(enaplot = g2.plot, points = g2.points.for.plot, colors = "red")
    sub.plot = ena.plot.points(enaplot = sub.plot, points = g1.points.for.plot, colors = "blue")
    sub.plot = ena.plot.points(enaplot = sub.plot, points = g2.points.for.plot, colors = "red")
  }

  if(mean == TRUE) {
    g1.points.for.plot = as.matrix(set$points)[group1.rows,,drop=FALSE]
    g2.points.for.plot = as.matrix(set$points)[group2.rows,,drop=FALSE]

    g1.plot = ena.plot.group(g1.plot, g1.points.for.plot, colors = "blue", labels = group1,confidence.interval = "box")
    g2.plot = ena.plot.group(g2.plot, g2.points.for.plot, colors = "red", labels = group2,confidence.interval = "box")
    sub.plot = ena.plot.group(sub.plot, g1.points.for.plot, colors = "blue", labels = group1,confidence.interval = "box")
    sub.plot = ena.plot.group(sub.plot, g2.points.for.plot, colors = "red", labels = group2,confidence.interval = "box")
  }

  else if(TRUE %in% c(network,points, mean) == FALSE) {
    stop("You must set at least one of points, mean, or network to TRUE to obtain a plot.")
  }

  set$plots[[group1]] = g1.plot
  set$plots[[group2]] = g2.plot
  set$plots[[paste0(group1,"-",group2)]] = sub.plot

  return(set)
}
"""