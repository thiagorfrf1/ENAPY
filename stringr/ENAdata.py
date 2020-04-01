ENAdata = """
####
#' ENAdata R6class
#'
#' @docType class
#' @importFrom R6 R6Class
#' @import data.table
#' @export
#'
#' @field raw A data frame constructed from the unit, convo, code, and metadata parameters of ena.accumulate.data
#' @field adjacency.vectors A data frame of adjacency (co-occurrence) vectors by row
#' @field accumulated.adjacency.vectors A data frame of adjacency (co-occurrence) vectors accumulated per unit
#' @field model The type of ENA model: EndPoint, Accumulated Trajectory, or Separate Trajectory
#' @field units A data frame of columns that were combined to make the unique units. Includes column for trajectory selections. (unique)
#' @field unit.names A vector of unique unit values
#' @field metadata A data frame of unique metadata for each unit
#' @field trajectories A list: units - data frame, for a given row tells which trajectory it's a part; step - data frame, where along the trajectory a row sits
#'
#' @field codes A vector of code names
#' @field function.call The string representation of function called and parameters provided
#' @field function.params A list of all parameters sent to function call
####
ENAdata <- R6::R6Class("ENAdata", public = list(
  initialize = function(
    file,
    units = NULL,
    units.used = NULL,
    units.by = NULL,
    conversations.by = NULL,
    codes = NULL,
    model = NULL,
    weight.by = "binary",
    window.size.back = 1,
    window.size.forward = 0,
    mask = NULL,
    include.meta = T,
    ...
  ) {
    args <- list(...);
    self$function.call <- sys.call(-1);
    self$function.params <- list();

    private$file <- file;
    self$units <- units;
    private$units.used <- units.used;
    private$units.by <- units.by
    private$conversations.by <- conversations.by;
    self$codes <- codes;

    if (is.data.frame(self$codes)) {
      self$codes <- colnames(self$codes);
    }

    private$weight.by <- weight.by;
    private$window.size <- list(
      "back" = window.size.back,
      "forward" = window.size.forward
    );

    for (p in c("units", "units.used", "units.by",
               "conversations.by", "codes", "model", "weight.by",
               "window.size.back", "window.size.forward", "mask",
               "in.par", "grainSize", "include.meta")
    ) {
      if (exists(x = p)) {
        self$function.params[[p]] <- get(p)
      } else if (!is.null(args[[p]])) {
        self$function.params[[p]] <- args[[p]]
      }
    }

    self$model <- model

    private$mask <- mask

    return(self)
  },

    ####
    ## Public Properties
    ####
      model = NULL,
      raw = NULL,
      adjacency.vectors = NULL,
      adjacency.matrix = NULL,
      accumulated.adjacency.vectors = NULL,
      adjacency.vectors.raw = NULL,
      units = NULL,
      unit.names = NULL,
      metadata = NULL,
      trajectories = list(
        units = NULL,
        step = NULL
      ),
      codes = NULL,
      function.call = NULL,
      function.params = NULL,
    ####
    ## END: Public Properties
    ####

    ####
    ## Public Functions
    ####
      ####
      # Process the accumulation
      #####
      process = function() {
        private$loadFile();
      },

      ####
      # \code{get()} - Return a read-only property
      # \preformatted{  Example:
      #     get( x = 'file' )}
      # \preformatted{  Parameters:
      #      x - Property to return. Defaults 'data', returns the original data
      ####
      get = function(x = "data") {
        return(private[[x]])
      },

      ####
      # \code{read()} - Return the accumulated data
      # \preformatted{  Example:
      #     get( colnames = T, sep = " & " )}
      # \preformatted{  Parameters:
      #      colnames - Logical, whether to replace colnames with their names
      #                 values from the adjacency (co-occurrence)
      #      sep - String to use as a seperator in the updated column names.
      #             Ignored if colnames == F}
      ####
      read = function(colnames = T, sep = " & ") {
        named_data <- data.table::copy(self$adjacency.vectors);
        if (colnames == T) {
          named_rows <- attr(named_data, "adjacency.matrix");
          colnames(named_data)[grep("adjacency.code", colnames(named_data))] <-
            apply(named_rows, 2, function(x) paste(x[1], x[2], sep = sep))
        }

        return()
      },
      add.metadata = function(merge = F) {
        meta_avail <- colnames(self$raw)[
          -which(colnames(self$raw) %in%
                  c(self$codes, private$units.by, private$conversations.by))]

        meta_avail <- meta_avail[which(meta_avail != "ENA_UNIT")]
        meta_cols_to_use <- meta_avail[apply(self$raw[, lapply(.SD, uniqueN),
                                                    by = c(private$units.by),
                                                    .SDcols = meta_avail
                                                 ][, c(meta_avail), with = F]
                                    , 2, function(x) all(x == 1))
                                  ]
        raw.meta <- self$raw[!duplicated(ENA_UNIT)][
                      ENA_UNIT %in% unique(
                        self$accumulated.adjacency.vectors$ENA_UNIT
                      ),
                      c("ENA_UNIT", private$units.by, meta_cols_to_use),
                      with = F
                    ]

        df_to_return <- NULL;
        if (merge == T) {
          df_to_return <- merge(
            self$adjacency.vectors,
            raw.meta[, unique(colnames(raw.meta)), with = F],
            by = c("ENA_UNIT"),
            suffixes = c("", ".y"), sort = F
          )
        } else {
          df_to_return <- raw.meta[ENA_UNIT %in% self$unit.names,];
        }

        return(df_to_return)
      },
      print = function(...) {
        args <- list(...);
        fields <- NULL;
        to.print <- list();
        if (is.null(args$fields)) {
          fields <- names(get(class(self))$public_fields)
        } else {
          fields <- args$fields
        }
        for (f in fields) {
          to.print[[f]] <- self[[f]]
        }
        return(to.print);
      }
    ####
    ## END: Public Functions
    ####
  ),

  ####
  ### Private
  ####
  private = list(

    ####
    ## Private Properties
    #####
      file = NULL,
      window.size = NULL,
      units.used = NULL,
      units.by = NULL,
      conversations.by = NULL,
      weight.by = NULL,
      trajectory.by = NULL,
      mask = NULL,
    #####
    ## END: Private Properties
    ####

    ####
    ## Private Functions
    #####
    loadFile = function() {
      if(any(class(private$file) == "data.table")) {
        df_DT <- private$file
      } else {
        if(any(class(private$file) == "data.frame")) {
          df <- private$file
        } else {
          df <- read.csv(private$file)
        }
        df_DT <- data.table::as.data.table(df)
      }

      self$raw <- data.table::copy(df_DT)
      self$raw$ENA_UNIT <- merge_columns_c(self$raw, private$units.by)

      self <- accumulate.data(self)
      self$units <- self$adjacency.vectors[, private$units.by, with = F]

      if (!self$model %in% c("AccumulatedTrajectory", "SeparateTrajectory")) {
        self$unit.names <- self$adjacency.vectors$ENA_UNIT
      } else {
        self$trajectories$units <- self$units
        conversation <- self$adjacency.vectors[,
                          private$conversations.by,
                          with = F
                        ]

        self$trajectories$step <- conversation
        self$units <- cbind(self$units, conversation)
        self$unit.names <- paste(
          self$adjacency.vectors$ENA_UNIT,
          self$adjacency.vectors$TRAJ_UNIT,
          sep = "."
        )
      }

      self$adjacency.vectors.raw <- self$adjacency.vectors

      adjCols <- colnames(self$adjacency.vectors)[
                  grep("adjacency.code", colnames(self$adjacency.vectors))
                ];

      if (is.null(private$mask)) {
        private$mask <- matrix(1,
                          nrow = length(self$codes),
                          ncol = length(self$codes),
                          dimnames = list(self$codes, self$codes))
      }

      self$adjacency.vectors[, c(adjCols)] <-
        self$adjacency.vectors[, c(adjCols), with = F] *
          rep(
            private$mask[upper.tri(private$mask)],
            rep(nrow(self$adjacency.vectors), length(adjCols))
          )

      if( is.function(private$weight.by) ) {
        cols <- colnames(self$adjacency.vectors)[
                  grep("adjacency.code", colnames(self$adjacency.vectors))
                ]
        self$adjacency.vectors <- self$adjacency.vectors[,
                                    lapply(
                                      .SD,
                                      private$weight.by
                                    ),
                                    .SDcols = cols,
                                    by = 1:nrow(self$adjacency.vectors)
                                  ]
      }

      if( self$function.params$include.meta == T) {
        self$metadata <- self$add.metadata(merge = F);
      } else {
        self$metadata <- data.frame();
      }

      self$adjacency.vectors <- self$adjacency.vectors[,
                                  grep("adjacency.code",
                                    colnames(self$adjacency.vectors)),
                                  with = F
                                ]

      return(self);
    }
  )
)
"""