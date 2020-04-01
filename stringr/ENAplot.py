ENAplot = """
####
#' ENAset R6class
#'
#' @docType class
#' @importFrom R6 R6Class
#' @import data.table
#' @export
#'
# @param enaset ENAplot Object
#'
#' @field enaset - The \code{\link{ENAset}} object from which the ENAplot was constructed
#' @field plot - The plotly object used for data visualization
#'
####

ENAplot = R6::R6Class("ENAplot",

  public = list(

    ####
    ## Constructor - documented in main class declaration
    ####
      initialize = function(
        enaset = NULL,

        title = "ENA Plot",

        dimension.labels = c("",""),

        font.size = 14,
        font.color = "#000000",
        font.family = "Arial",
        scale.to = "network",
        ...
      ) {

        if (is(enaset, "ENAset")) {
          warning(paste0("Usage of ENAset objects will be deprecated ",
            "and potentially removed altogether in future versions."))

          enaset <- ena.set(enaset);
        }

        code.cols = !colnames(enaset$line.weights) %in% colnames(enaset$meta.data)

        args = list(...);
        if(!is.null(args$multiplier)) {
          private$multiplier = args$multiplier
        }
        if(!is.null(args$point.size)) {
          self$point$size = args$point.size
        }
        self$enaset <- enaset;

        private$title <- title;
        private$dimension.labels <- dimension.labels;
        private$font.size <- font.size;
        private$font.color <- font.color;
        private$font.family <- font.family;
        private$font = list (
          size = private$font.size,
          color = private$font.color,
          family = private$font.family
        );
        self$plot <- plotly::plot_ly(
          mode = "markers",
          type ="scatter"
        );

        if (is.list(scale.to)) {
          max.axis = max(abs(as.matrix(enaset$points)))*1.2
          if(is.null(scale.to$x)) {
            axis.range.x = c(-max.axis, max.axis)
          } else {
            axis.range.x = scale.to$x
          }
          if(is.null(scale.to$y)) {
            axis.range.y = c(-max.axis, max.axis)
          } else {
            axis.range.y = scale.to$y
          }
        } else {
          if(scale.to == "points") {
            max.axis = max(abs(as.matrix(enaset$points)))*1.2
          } else if (is.numeric(scale.to)) {
            max.axis = tail(scale.to, 1)
          } else {
            max.axis = max(abs(as.matrix(enaset$rotation$nodes)))*1.2;
          }
          axis.range.x = axis.range.y = c(-max.axis, max.axis)
        }

        graph.axis <- list(
          titlefont = private$font,
          showgrid = F,
          zeroline = T,
          showticklabels = T,
          showgrid = T
          # range=c(-max.axis,max.axis)
        );
        if(!is.null(args$ticks)) {
          graph.axis$showticklabels = T;
          graph.axis$ticks = args$ticks$location;
          graph.axis$tickcolor = args$ticks$color;
          graph.axis$tickangle = args$ticks$angle;
        }
        self$axes$x = graph.axis
        self$axes$x$title = dimension.labels[1];
        self$axes$x$range = axis.range.x
        self$axes$y = graph.axis
        self$axes$y$title = dimension.labels[2];
        self$axes$y$range = axis.range.y

        self$plot = plotly::layout(
          self$plot,
          title =  title,
          xaxis = self$axes$x,
          yaxis = self$axes$y,
          font = list (
            size = 12,
            color = private$font.color,
            family = private$font.family
          )
        );
      },
    ####
    ## END: Constructor - documented in main class declaration
    ####

    ####
    ## Public Properties
    ####
      enaset = NULL,
      plot = NULL,
      axes = list(
        x = NULL, y = NULL
      ),
      point = list(
        size = 5
      ),
      palette = c("#386CB0", "#F0027F", "#7FC97F", "#BEAED4",
                  "#FDC086","#FFFF99", "#BF5B17"),
      plotted = list(
        points = list(), networks = list(),
        trajectories = list(), means = list()
      ),
    ####
    ## END: Public Properties
    ####

    ####
    ## Public Functions
    ####
      print = function() {
        print(self$plot);
      },

      ####
      # \code{get()} - Return a read-only property
      # \preformatted{  Example:
      #     get( x = 'title' )}
      # \preformatted{  Parameters:
      #      x - Property to return. Defaults to 'title', returning the title}
      ####
      get = function(x) {
        return(private[[x]])
      },
      update.axis.lines = function() {

      }
    ####
    ## END: Public Functions
    ####
  ),

  private = list(
    ####
    ## Private Properties
    ####
      title = "ENA Plot",

      dimension.labels = c("X","Y"),

      font = list(),
      font.size = 14,
      font.color = "#000000",
      font.family = "Arial",
      #plot.color = I("black"),

      multiplier = 5
    ####
    ## END: Private Properties
    ####
  )
)
"""