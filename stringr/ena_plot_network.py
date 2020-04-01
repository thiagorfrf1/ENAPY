ena_plot_network = """
##
#' @title Plot an ENA network
#'
#' @description Plot an ENA network: nodes and edges
#'
#' @details lots a network graph, including nodes (taken from codes in the ENAplot) and the edges (provided in network)
#'
#' @export
#'
#' @param enaplot \code{\link{ENAplot}} object to use for plotting
#' @param network dataframe or matrix containing the edge weights for the network graph; typically comes from ENAset$line.weights
#' @param node.positions matrix containing the positiions of the nodes. Defaults to enaplot$enaset$node.positions
#' @param adjacency.key matrix containing the adjacency key for looking up the names and positions
#' @param colors A String or vector of colors for positive and negative line weights. E.g. red or c(pos= red, neg = blue), default: c(pos= red, neg = blue)
#' @param edge_type A String representing the type of line to draw, either "line", "dash", or "dot"
#' @param show.all.nodes A Logical variable, default: true
#' @param threshold A vector of numeric min/max values, default: c(0,Inf) plotting . Edge weights below the min value will not be displayed; edge weights above the max value will be shown at the max value.
#' @param thin.lines.in.front A logical, default: true
#' @param thickness A vector of numeric min/max values for thickness, default:  c(min(abs(network)), max(abs(network)))
#' @param opacity A vector of numeric min/max values for opacity, default: thickness
#' @param saturation A vector of numeric min/max values for saturation, default: thickness
#' @param scale.range A vector of numeric min/max to scale from, default: c(0.1,1) or if min(network) is 0, c(0,1)
#' @param node.size A lower and upper bound used for scaling the size of the nodes, default c(0, 20)
#' @param labels A character vector of node labels, default: code names
#' @param label.offset A character vector of representing the positional offset relative to the respective node. Defaults to "middle right" for all nodes. If a single values is provided, it is used for all positions, else the length of the
#' @param label.font.size An integer which determines the font size for graph labels, default: enaplot$font.size
#' @param label.font.color A character which determines the color of label font, default: enaplot$font.color
#' @param label.font.family A character which determines font type, choices: Arial, Courier New, Times New Roman, default: enaplot$font.family
#' @param legend.name A character name used in the plot legend. Not included in legend when NULL (Default), if legend.include.edges is TRUE will always be "Nodes"
#' @param legend.include.edges Logical value indicating if the edge names should be included in the plot legend. Forces legend.name to be "Nodes"
#' @param scale.weights Logical indicating to scale the supplied network
#' @param ... Additional parameters
#'
#' @seealso \code{\link{ena.plot}}, \code{\link{ena.plot.points}}
#' @importFrom scales rescale

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
#'   enadata = accum,
#'   rotation.by = ena.rotate.by.mean,
#'   rotation.params = list(
#'     accum$meta.data$Condition=="FirstGame",
#'     accum$meta.data$Condition=="SecondGame"
#'   )
#' )
#'
#' plot = ena.plot(set)
#'
#' ### Subset rotated points and plot Condition 1 Group Mean
#' as.matrix(set$points$Condition$FirstGame)
#'
#' first.game.points = as.matrix(set$points$Condition$FirstGame)
#' plot = ena.plot.group(plot, first.game.points, labels = "FirstGame",
#'     colors = "red", confidence.interval = "box")
#'
#' ### Subset rotated points and plot Condition 2 Group Mean
#' second.game.points = as.matrix(set$points$Condition$SecondGame)
#' plot = ena.plot.group(plot, second.game.points, labels = "SecondGame",
#'     colors  = "blue", confidence.interval = "box")
#'
#' ### get mean network plots
#' first.game.lineweights = as.matrix(set$line.weights$Condition$FirstGame)
#' first.game.mean = colMeans(first.game.lineweights)
#'
#' second.game.lineweights = as.matrix(set$line.weights$Condition$SecondGame)
#' second.game.mean = colMeans(second.game.lineweights)
#'
#' subtracted.network = first.game.mean - second.game.mean
#' plot = ena.plot.network(plot, network = subtracted.network)
#' print(plot)
#'
#' @return The \code{\link{ENAplot}} provided to the function, with its plot updated to include the nodes and provided connecting lines.
##
ena.plot.network = function(
  enaplot = NULL,
  network = NULL,
  node.positions = enaplot$enaset$rotation$nodes,
  adjacency.key = NULL, #enaplot$enaset$enadata$adjacency.matrix,
  colors = c(pos=enaplot$palette[1], enaplot$palette[2]),
  edge_type = "line", #c("line", "dash", "dot"),
  show.all.nodes = T,
  threshold = c(0),
  thin.lines.in.front = T,

  thickness = c(min(abs(network)), max(abs(network))),
  opacity = thickness,
  saturation = thickness,
  scale.range = c(ifelse(min(network)==0, 0, 0.1), 1),

  node.size = c(3,10),

  labels = NULL,
  label.offset = "middle right",
  label.font.size = enaplot$get("font.size"),
  label.font.color = enaplot$get("font.color"),
  label.font.family = enaplot$get("font.family"),
  legend.name = NULL,
  legend.include.edges = F,
  scale.weights = T,
  ...
) {
  if(choose(nrow(node.positions), 2) != length(network)) {
    stop(paste0("Network vector needs to be of length ", choose(nrow(node.positions), 2)))
  }
  node.rows <- NULL
  if(is(node.positions, "ena.nodes")) {
    if(is.null(adjacency.key)) {
      adjacency.key <- namesToAdjacencyKey(node.positions$code)
    }
    node.rows <- node.positions$code

    if(is.null(labels)) {
      labels <- node.positions$code
    }
  } else {
    if(is.matrix(node.positions)) {
      node.positions <- as.data.frame(node.positions)
    }
    adjacency.key <- namesToAdjacencyKey(rownames(node.positions))
    node.rows <- rownames(node.positions)
    if(is.null(labels)) {
      labels  <- rownames(node.positions)
    }
  }
  args = list(...);
  network.edges.shapes = list();
  edge_type = match.arg(arg = edge_type, choices = c("line", "dash", "dot"));

  nodes = data.frame(as.matrix(node.positions));
  colnames(nodes) = paste0("X", seq(colnames(nodes)))
  nodes$weight = rep(0, nrow(nodes))
  nodes$color = "black";

  # Handle label parameters
  if(length(label.offset) == 1) {
    label.offset = rep(label.offset[1], length(labels))
  }
  if(length(label.offset) != length(labels)) {
    stop("length(label.offset) must be equal to 1 or length(labels)")
  }

  # Handle legend parameters
  if(legend.include.edges == T && !is.null(legend.name)) {
    legend.name = "Nodes"
  }

  network.scaled = network;
  if(!is.null(threshold)) {
    multiplier.mask = ((network.scaled >= 0) * 1) - ((network.scaled < 0) * 1)
    if(length(threshold) == 1) {
      threshold[2] = Inf;
    } else if(threshold[2] < threshold[1]) {
      stop("Minimum threshold value must be less than the maximum value.");
    }

    if(threshold[1] > 0) {
      # network.scaled = network.scaled[sizes > threshold[1]]
      network.scaled[abs(network.scaled) < threshold[1]] = 0
    }
    if(threshold[2] < Inf && any(abs(network.scaled) > threshold[2]))  {
      to.threshold = abs(network.scaled) > threshold[2]
      network.scaled[to.threshold] = threshold[2]
      network.scaled[to.threshold] = network.scaled[to.threshold] * multiplier.mask[to.threshold]
    }
  }
  network.thickness = abs(network.scaled);
  network.saturation = abs(network.scaled);
  network.opacity = abs(network.scaled);

  network.to.keep = (network != 0) * 1
  if(!is.null(args$scale.weights) && args$scale.weights == T) {
    network.scaled = network * (1 / max(abs(network)));

    network.thickness = scales::rescale(x = abs(network.scaled), to = scale.range, from = thickness);
  }
  network.scaled = network.scaled * network.to.keep
  network.thickness = network.thickness * network.to.keep

  network.saturation = scales::rescale(x = abs(network.scaled), to = scale.range, from = saturation);
  network.opacity = scales::rescale(x = abs(network.scaled), to = scale.range, from = opacity);

  pos.inds = as.numeric(which(network.scaled >=0));
  neg.inds = as.numeric(which(network.scaled < 0));

  colors.hsv = rgb2hsv(col2rgb(colors))

  if(ncol(colors.hsv) == 1) {
    colors.hsv[[4]] = colors.hsv[1] + 0.5;
    if(colors.hsv[4] > 1) {
      colors.hsv[4] = colors.hsv[4] - 1;
    }

    colors.hsv[[5]] = colors.hsv[2];
    colors.hsv[[6]] = colors.hsv[3];
    dim(colors.hsv) = c(3,2);
  }

  mat = as.matrix(adjacency.key);
  for (i in 1:length(network)) {
    v0 <- nodes[node.rows==mat[1,i], ];
    v1 <- nodes[node.rows==mat[2,i], ];
    nodes[node.rows==mat[1,i],]$weight = nodes[node.rows==mat[1,i],]$weight + abs(network.thickness[i]);
    nodes[node.rows==mat[2,i],]$weight = nodes[node.rows==mat[2,i],]$weight + abs(network.thickness[i]);

    color = NULL
    if(i %in% pos.inds) {
      color = colors.hsv[,1];
    } else {
      color = colors.hsv[,2];
    }
    color[2] = network.saturation[i];

    edge_shape = list(
      type = "line",
      opacity = network.opacity[i],
      nodes = c(mat[,i]),
      line = list(
        name = "test",
        color= hsv(color[1],color[2],color[3]),
        width= abs(network.thickness[i]) * enaplot$get("multiplier"),
        dash = edge_type
      ),
      x0 = as.numeric(v0[1]),
      y0 = as.numeric(v0[2]),
      x1 = as.numeric(v1[1]),
      y1 = as.numeric(v1[2]),
      layer = "below",
      size = as.numeric(abs(network.scaled[i]))
    );
    network.edges.shapes[[i]] = edge_shape
  };

  if(thin.lines.in.front) {
    network.edges.shapes = network.edges.shapes[rev(order(sapply(network.edges.shapes, "[[", "size")))]
  } else {
    network.edges.shapes = network.edges.shapes[order(sapply(network.edges.shapes, "[[", "size"))]
  }

  rows.to.keep = rep(T, nrow(nodes))
  if(show.all.nodes == F) {
    rows.to.keep = nodes$weight != 0
    # nodes = nodes[rownames(nodes) %in% unique(as.character(sapply(network.edges.shapes, "[[", "nodes"))), ]
  }
  nodes = nodes[rows.to.keep,];
  mode = "markers+text"
  if(!is.null(args$labels.hide) && args$labels.hide == T) {
    mode="markers"
  }
  nodes$weight = scales::rescale((nodes$weight * (1 / max(abs(nodes$weight)))), node.size) # * enaplot$get("multiplier"));

  show.legend = !is.null(legend.name);
  if(legend.include.edges) {
    if(is.null(legend.name)) {
      legend.name = "Nodes"
    }
    show.legend = T;
  }

  enaplot$plot = plotly::add_trace(
    enaplot$plot,
    type = "scatter",
    data = nodes,
    x = ~X1,
    y = ~X2,
    mode = mode,
    textposition = label.offset[rows.to.keep],
    marker = list(
      color = "#000000",
      size = abs(nodes$weight)
      #,name = labels[i] #rownames(nodes)[i]
    ),
    textfont = list (
      family = label.font.family,
      size = label.font.size,
      color = label.font.color
    ),
    text = labels[rows.to.keep], #rownames(nodes),
    legendgroup = legend.name,
    name = legend.name,
    showlegend = show.legend,
    hoverinfo = 'none'
  );

  if (length(network.edges.shapes) > 0 ) {
    enaplot$plotted$networks[[length(enaplot$plotted$networks) + 1]] <- network.edges.shapes

    for (n in 1:length(network.edges.shapes)) {
      e = network.edges.shapes[[n]];

      name = NULL;
      show.legend = F;
      this.name = paste(e$nodes[1],e$nodes[2], sep=".")
      if(legend.include.edges) {
        name = this.name;
        show.legend = T;
      }

      enaplot$plot = plotly::add_trace(
        enaplot$plot,
        type = "scatter",
        mode = "lines",
        data = data.frame(X1=c(e$x0,e$x1), X2=c(e$y0,e$y1)),
        x = ~X1, y = ~X2,
        line = e$line,
        opacity = e$opacity,
        legendgroup = if(legend.include.edges == T) this.name else legend.name,
        showlegend = show.legend,
        name = name
      )
    }
  }

  enaplot
}
"""