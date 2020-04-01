ena_plot_points = """
##
#' @title Plot points on an ENAplot
#'
#' @description Plot all or a subset of the points of an ENAplot using the plotly plotting library
#'
#' @export
#'
#' @param enaplot \code{\link{ENAplot}} object to use for plotting
#' @param points A dataframe of matrix where the first two column are X and Y coordinates
#' @param point.size A data.frame or matrix where the first two column are X and Y coordinates of points to plot in a projected ENA space defined in ENAplot
#' @param labels A character vector of point labels, length nrow(points); default: NULL
#' @param confidence.interval A character determining markings to use for confidence intervals, choices: none, box, crosshair, default: none
#' @param outlier.interval A character determining markings to use for outlier interval, choices: none, box, crosshair, default: none
#' @param confidence.interval.values A matrix/dataframe where columns are CI x and y values for each point
#' @param outlier.interval.values A matrix/dataframe where columns are OI x and y values for each point
#' @param shape A character which determines the shape of point markers, choices:   square, triangle, diamond, circle, default: circle
#' @param colors A character vector of the point marker colors; if one given it is used for all, otherwise must be same length as points; default: black
#' @param label.offset character: top left (default), top center, top right, middle left, middle center, middle right, bottom left, bottom center, bottom right
#' @param label.group A string used to group the labels in the legend. Items plotted with the same label.group will show/hide together when clicked within the legend.
#' @param label.font.size An integer which determines the font size for point labels, default: enaplot$font.size
#' @param label.font.color A character which determines the color of label font, default: enaplot$font.color
#' @param label.font.family	A character which determines label font type, choices: Arial, Courier New, Times New Roman, default: enaplot$font.family
#' @param show.legend Logical indicating whether to show the point labels in the in legend
#' @param legend.name Character indicating the name to show above the plot legend
#' @param texts [TBD]
#' @param ... additional parameters addressed in inner function
#'
#'
#' @seealso \code{\link{ena.plot}}, \code{\link{ENAplot}}, \code{\link{ena.plot.group}}
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
#'   enadata = accum,
#'   rotation.by = ena.rotate.by.mean,
#'   rotation.params = list(
#'       accum$meta.data$Condition=="FirstGame",
#'       accum$meta.data$Condition=="SecondGame"
#'   )
#' )
#'
#' plot = ena.plot(set)
#'
#' group1.points = set$points[set$meta.data$Condition == "FirstGame",]
#' group2.points = set$points[set$meta.data$Condition == "SecondGame",]
#' plot = ena.plot.points(plot, points = group1.points);
#' plot = ena.plot.points(plot, points = group2.points);
#' print(plot);
#'
#' @return \code{\link{ENAplot}} The ENAplot provided to the function, with its plot updated to include the new points.
##
ena.plot.points = function(
  enaplot,

  points = NULL,    #vector of unit names or row indices
  point.size = enaplot$point$size,
  labels = NULL, #unique(enaplot$enaset$enadata$unit.names),
  label.offset = "top left",
  label.group = NULL,

  label.font.size = NULL, #enaplot$get("font.size"),
  label.font.color = NULL, #enaplot$get("font.color"),
  label.font.family = NULL, #enaplot$get("font.family"),

  shape = "circle",
  colors = NULL, # c("blue"), #rep(I("black"), nrow(points)),

  confidence.interval.values = NULL,
  confidence.interval = c("none", "crosshairs", "box"),

  outlier.interval.values = NULL,
  outlier.interval = c("none", "crosshairs", "box"),
  show.legend = T,
  legend.name = "Points",
  texts = NULL,
  ...
) {
  ###
  # Parameter Checking and Cleaning
  ###
    env = environment();
    for(n in c("font.size", "font.color", "font.family")) {
      if(is.null(get(paste0("label.",n))))
        env[[paste0("label.",n)]] = enaplot$get(n);
    }

    if(is.null(points)) {
      # stop("Must provide points to plot.")
      points = enaplot$enaset$points
    }

    if(is(points, "numeric")){
      points = matrix(points);
      dim(points) = c(1,nrow(points))
      points.layout = data.table::data.table(points);
    }
    else if (is.data.table(points)) {
      # points.layout = remove_meta_data(points)
      points.layout = data.table::copy(points)
    }
    else {
      points.layout = data.table::data.table(points);
    }

    if(!is.character(label.font.family)) {
      label.font.family = enaplot$get("font.family");
    }

    confidence.interval = match.arg(confidence.interval);
    outlier.interval = match.arg(outlier.interval);

    # shape = match.arg(shape);
    valid.shapes = c("circle", "square", "triangle-up", "diamond");
    if(!all(shape %in% valid.shapes))
      stop(sprintf( "Unrecognized shapes: %s", paste(unique(shape[!(shape %in% valid.shapes)]), collapse = ", ") ))
    if(length(shape) == 1)
      shape = rep(shape, nrow(points.layout))

    valid.label.offsets = c("top left","top center","top right","middle left","middle center","middle right","bottom left","bottom center","bottom right");
    if(!all(label.offset %in% valid.label.offsets))
      stop(sprintf( "Unrecognized label.offsets: %s", paste(unique(label.offset[!(label.offset %in% valid.label.offsets)]), collapse = ", ") ))
    if(length(label.offset) == 1)
      label.offset = rep(label.offset, nrow(points.layout))

    if(grepl("^c", confidence.interval) && grepl("^c", outlier.interval)) {
      print("Confidence Interval and Outlier Interval cannot both be crosshair");
      print("Plotting Outlier Interval as box");
      outlier.interval = "box";
    }

    if(length(colors) == 1) {
      colors = rep(colors, nrow(points.layout))
    }
    if(length(point.size) == 1)
      point.size = rep(point.size, nrow(points.layout))
    if(is.null(labels))
      show.legend = F
  ###
  # END: Parameter Checking and Cleaning
  ###

  ###
  # Set error value for CI|OI crosshair on plot
  ###
    error = list(x = list(visible=T, type="data"), y = list(visible=T, type="data"));
    int.values = NULL;
    if(grepl("^c", confidence.interval) && !is.null(confidence.interval.values)) {
      int.values = confidence.interval.values;
    }
    else if(grepl("^c", outlier.interval) && !is.null(outlier.interval.values)) {
      int.values = outlier.interval.values;
    }
    error$x$array = int.values[,1];
    error$y$array = int.values[,2];
  ###
  # END: Set error value for crosshair on plot
  ###

  ###
  # Set box value for CI|OI box on plot
  #####
    box.values = NULL;
    if(grepl("^b", confidence.interval) && !is.null(confidence.interval.values)) {
      box.values = confidence.interval.values;
      box.label = "Conf. Int.";
    }
    if(grepl("^b", outlier.interval) && !is.null(outlier.interval.values)) {
      box.values = outlier.interval.values;
      box.label = "Outlier Int.";
    }
  ######
  # END: Set box value for CI|OI box on plot
  ###

  ###
  # Plot
  #####
    points.matrix = remove_meta_data(points.layout)
    colnames(points.matrix) = paste0("X", rep(1:ncol(points.matrix)));
    this.max = max(points.matrix);
    for(m in 1:nrow(points.matrix)) {
      enaplot$plot = plotly::add_trace(
        p = enaplot$plot,
        data = points.matrix[m,],
        type ="scatter",
        x = ~X1, y = ~X2,
        mode = "markers+text",
        marker = list(
          symbol = shape[m],
          color = colors[m],
          size = point.size[m]
        ),
        error_x = error$x, error_y = error$y,
        showlegend = show.legend,
        # legendgroup = label.group,
        # legendgroup = ifelse(!is.null(box.label), labels[1], NULL),
        name = labels[m],
        text = texts[m], #labels[m],
        textfont = list(
          family = label.font.family,
          size = label.font.size,
          color = label.font.color
        ),
        legendgroup = legend.name,
        textposition = label.offset[m],
        hoverinfo = "x+y+name"
      )
    }

    if(!is.null(box.values)) {
      boxv = data.frame(
        X1 = c(box.values[1,1], box.values[2,1], box.values[2,1], box.values[1,1] ,box.values[1,1]),
        X2 = c(box.values[1,2], box.values[1,2], box.values[2,2], box.values[2,2], box.values[1,2])
      )
      this.max = max(boxv, this.max)
      enaplot$plot = plotly::add_trace(
        p = enaplot$plot,
        data = boxv,
        type = "scatter",
        x = ~X1, y = ~X2,
        mode = "lines",
        line = list(
          width = 1,
          color = colors[1],
          dash = "dash"
        ),
        # "legendgroup" = labels[1],
        showlegend = show.legend,
        name = box.label
      )
    }

    if(this.max*1.2 > max(enaplot$axes$y$range)) {
      this.max = this.max * 1.2
      enaplot$axes$x$range = c(-this.max, this.max)
      enaplot$axes$y$range = c(-this.max, this.max)
      enaplot$plot = plotly::layout(
        enaplot$plot,
        xaxis = enaplot$axes$x,
        yaxis = enaplot$axes$y
      );
    }
  #####
  # END: Plot
  ###

  return(enaplot);
}

"""