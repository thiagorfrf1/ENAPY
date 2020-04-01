ena_plot_trajectory = """
#####
#' @title Plot of ENA trajectories
#'
#' @description Function used to plot trajectories
#'
#' @export
#'
#' @param enaplot \code{\link{ENAplot}} object to use for plotting
#' @param points dataframe of matrix - first two column are X and Y coordinates, each row is a point in a trajectory
#' @param by vector used to subset points into individual trajectories, length nrow(points)
#' @param names character vector - labels for each trajectory of points, length length(unique(by))
#' @param labels character vector - point labels, length nrow(points)
#' @param labels.show A character choice: Always, Hover, Both.  Default: Both
# @param confidence.interval A character that determines which confidence interval type to use, choices: none, box, crosshair, default: none
# @param outlier.interval A character that determines which outlier interval type to use, choices: none, box, crosshair, default: none
# @param confidence.interval.values A matrix/dataframe where columns are CI x and y values for each point
# @param outlier.interval.values A matrix/dataframe where columns are OI x and y values for each point
#' @param colors A character vector, that determines marker color, default NULL results in
#' alternating random colors. If single color is supplied, it will be used for all
#' trajectories, otherwise the length of the supplied color vector should be equal
#' to the length of the supplied names (i.e a color for each trajectory being plotted)
#' @param shape A character which determines the shape of markers, choices: square, triangle, diamond, circle, default: circle
#' @param label.offset A numeric vector of an x and y value to offset labels from the coordinates of the points
#' @param label.font.size An integer which determines the font size for labels, default: enaplot\$font.size
#' @param label.font.color A character which determines the color of label font, default: enaplot\$font.color
#' @param label.font.family A character which determines font type, choices: Arial, Courier New, Times New Roman, default: enaplot\$font.family
#' @param default.hidden A logical indicating if the trajectories should start hidden (click on the legend to show them) Default: FALSE
#'
#' @seealso \code{\link{ena.plot}}
#'
#' @examples
#' data(RS.data)
#'
#' codeNames = c('Data','Technical.Constraints','Performance.Parameters',
#'   'Client.and.Consultant.Requests','Design.Reasoning','Collaboration');
#'
#' accum = ena.accumulate.data(
#'   units = RS.data[,c("UserName","Condition")],
#'   conversation = RS.data[,c("GroupName","ActivityNumber")],
#'   metadata = RS.data[,c("CONFIDENCE.Change","CONFIDENCE.Pre","CONFIDENCE.Post","C.Change")],
#'   codes = RS.data[,codeNames],
#'   window.size.back = 4,
#'   model = "A"
#' );
#'
#' set = ena.make.set(accum);
#'
#' ### get mean network plots
#' first.game.lineweights = as.matrix(set$line.weights$Condition$FirstGame)
#' first.game.mean = colMeans(first.game.lineweights)
#'
#' second.game.lineweights = as.matrix(set$line.weights$Condition$SecondGame)
#' second.game.mean = colMeans(second.game.lineweights)
#'
#' subtracted.network = first.game.mean - second.game.mean
#'
#' # Plot dimension 1 against ActivityNumber metadata
#' dim.by.activity = cbind(
#'     as.matrix(set$points)[,1],
#'     set$trajectories$ActivityNumber * .8/14-.4  #scale down to dimension 1
#' )
#'
#' plot = ena.plot(set)
#' plot = ena.plot.network(plot, network = subtracted.network, legend.name="Network")
#' plot = ena.plot.trajectory(
#'   plot,
#'   points = dim.by.activity,
#'   names = unique(set$model$unit.label),
#'   by = set$trajectories$ENA_UNIT
#' );
#' print(plot)
#'
#' @return The \code{\link{ENAplot}} provided to the function, with its plot updated to include the trajectories
#####
ena.plot.trajectory = function(
  enaplot,
  points,
  by = NULL,
  labels = NULL, #unique(enaplot$enaset$enadata$units),
  labels.show = c("Always","Hover","Both"),
  names = NULL,
  label.offset = NULL,
  label.font.size = enaplot$get("font.size"),
  label.font.color = enaplot$get("font.color"),
  label.font.family = c("Arial", "Courier New", "Times New Roman"),
  shape = c("circle", "square", "triangle-up", "diamond"),
  colors = NULL,
  default.hidden = F
) {
  if(!is.character(label.font.family)) {
    label.font.size = enaplot$get("font.family");
  }
  labels.show <- match.arg(labels.show);
  shape <- match.arg(shape);

  if(is.null(by)) {
    by <- list(all = rep(T, nrow(points)));
  }
  if(!is(points, "data.table")) {
    points <- data.table::as.data.table(points);
  }
  if(length(colors) == 1)
    colors <- rep(colors, length(names))

  mode <- "lines+markers+text";
  hoverinfo <- "x+y";
  tbl <- data.table::data.table(points);
  if (!is.null(labels)) {
    if (labels.show %in% c("Always","Both"))
      mode <- paste0(mode,"+text");
    if (labels.show %in% c("Hover","Both"))
      hoverinfo <- paste0(hoverinfo,"+text");

    tbl = data.table::data.table(points, labels = labels);
  }

  if(!is.null(by)) {
    if(is.character(by) && length(by) == nrow(tbl))
        by <- as.factor(by)

    dfdt_trajs <- tbl[,{ data.table::data.table(lines = list(.SD))  }, by = by]
  } else {
    dfdt_trajs <- tbl[,{ data.table::data.table(lines = list(.SD))  }]
  }

  valid_label_offsets = c("top left","top center","top right","middle left",
              "middle center","middle right","bottom left","bottom center",
              "bottom right")
  if(!all(label.offset %in% valid_label_offsets))
    stop(sprintf( "Unrecognized label.offsets: %s",
      paste(unique(label.offset[!(label.offset %in% valid_label_offsets)]),
      collapse = ", ") ))

  if(length(label.offset) == 1)
    label.offset = rep(label.offset, nrow(dfdt_trajs))

  if (!is.null(colors) &&
      length(colors) > 1 && length(colors) != length(names)
  ) {
    stop("Length of the colors must be 1 or the same length as by")
  }

  for (x in 1:nrow(dfdt_trajs)) {
    d <- remove_meta_data(dfdt_trajs[x,]$lines[[1]])
    d.names <- colnames(d)
    enaplot$plot = plotly::add_trace(
      enaplot$plot,
      data = d,
      x = as.formula(paste0("~", d.names[1])),
      y = as.formula(paste0("~", d.names[2])),
      name = names[x],
      mode = mode,
      text = dfdt_trajs[x,]$lines[[1]]$labels,
      textposition = label.offset[x],
      hoverinfo = hoverinfo,
      showlegend = T,
      line = list (
        color = if(!is.null(colors)) colors[x] else NULL
      ),
      marker = list (
        symbol = shape
        ,color = if(!is.null(colors)) colors[x] else NULL
      ),
      textfont = list (
        family = label.font.family,
        size = label.font.size,
        color = label.font.color
      ),
      visible = ifelse(default.hidden, "legendonly", T)
    );
  }

  enaplot$plotted$trajectories[[
    length(enaplot$plotted$trajectories) + 1
  ]] <- dfdt_trajs

  return(enaplot);
}
"""