# -*- coding: UTF-8 -*-
import rpy2.robjects as robjects

from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import SignatureTranslatedAnonymousPackage as STAP
from rpy2.robjects import IntVector, Formula
import pandas as pd
pandas2ri.activate()



accumulate_data = """
accumulate.data <- function(enadata) {
  dfDT <- enadata$raw;

  units.used <- enadata$get("units.used")
  units.by <- enadata$get("units.by")
  trajectory.by <- enadata$get("trajectory.by")
  codes <- enadata$codes

  if (is.data.frame(codes)) {
    codes <- colnames(codes);
  }

  conversations.by <- enadata$get("conversations.by")
  window <- enadata$get("window.size")
  binaryStanzas <- F
  units.exclude <- enadata$get("units.exclude")

  if(is.null(trajectory.by)) {
    trajectory.by = conversations.by
  }

  ### should work to determine if binary is desired
  binary <- T;
  if (!identical(enadata$get("weight.by"), "binary")) {
    binary <- F
  } else {
    binary <- T
  }

  ### We need data
  if (is.null(dfDT) || nrow(dfDT) < 1) {
    stop("The provided data is NULL")
  }

  ###
  # We need a data.table, it's worth it.
  ###
  if(!data.table::is.data.table(dfDT)) {
    dfDT <- data.table::as.data.table(dfDT)
  }

  ###
  # Make a copy of the data for safe usage
  ###
  dfDT_codes <- data.table::copy(dfDT)

  ###
  # Create a column representing the ENA_UNIT as defined
  # by the the `units.by` parameter
  ###
  if(!"ENA_UNIT" %in% colnames(dfDT_codes)) {
    dfDT_codes$ENA_UNIT <- enadata$raw$ENA_UNIT <- merge_columns_c(
      dfDT_codes,
      cols = units.by, sep = "."
    )
  }

  ##
  # String vector of codesnames representing the names of the co-occurrences
  ##
  vL <- length(codes);
  adjacency.length <- ( (vL * (vL + 1)) / 2) - vL ;
  codedTriNames <- paste("adjacency.code",rep(1:adjacency.length), sep=".");

  initial_cols <- c(units.by, codes)
  just_codes <- c(codes)

  ##
  # Accumulated windows appended to the end of each row
  #
  # FIXME: Don't append on the results to the initial data.table,
  #        keep a separate to lookup the results for the co-occurred
  #        values later on.
  ##
  if (window$back == 1 && window$forward == 0) {
    dfDT.co.occurrences <- dfDT_codes[,{
        ocs <- data.table::as.data.table(
                rows_to_co_occurrences(
                  .SD[,.SD,.SDcols=codes, with=T],
                  binary = binary
                )
              );

        # Return value from data.table back to dfDT.co.occurrences
        data.table::data.table(.SD, ocs)
      },
      .SDcols = c(codes, conversations.by, trajectory.by),
      with = T
    ]

    ### Generate the ENA_UNIT column
    dfDT.co.occurrences$ENA_UNIT <- dfDT_codes$ENA_UNIT

    ### Keep original columns used for units
    dfDT.co.occurrences[, (units.by) := dfDT_codes[, .SD, .SDcols = units.by]]
  }
  else if (window$back == "Conversation") {
    ###
    # First sum all lines by conversation and unit to get vectors of codes
    # occurring in the whole conversation for each unit
    ###
    dfDT.conv.sum <- dfDT_codes[,
      lapply(.SD, sum), by = c(unique(conversations.by)),
      .SDcols = c(codes),
      with = T
    ]

    ###
    # Convert each units converstation sums into adjacency vectors
    ###
    # browser()
    dfDT.co.occurrences <- dfDT.conv.sum[,{
        ocs = data.table::as.data.table(rows_to_co_occurrences(.SD[,.SD,.SDcols=codes, with=T], binary = binary));
        data.table::data.table(.SD,ocs, ENA_UNIT=merge_columns_c(.SD, cols = units.by, sep="."))
      },
      .SDcols=unique(c(codes, conversations.by, trajectory.by, units.by)),
      with=T
    ];
  }
  else {
    ## parallell: https://stackoverflow.com/questions/14759905/data-table-and-parallel-computing
    ### Calculate occurrences of code within the provided window

    # if(enadata$function.params$in.par == T) {
    #   grainSize = ifelse(!is.null(enadata$function.params$grainSize), enadata$function.params$grainSize, 10);
    #   dfDT.co.occurrences = dfDT_codes[,
    #                            (codedTriNames) := try_one(
    #                              .SD[,.SD, .SDcols=just_codes],
    #                              window=window$back,
    #                              binary = binary,
    #                              grainSize = grainSize
    #                            ),
    #                            by=conversations.by,
    #                            .SDcols=initial_cols,
    #                            with=T
    #                          ];
    #
    # } else {
      dfDT.co.occurrences <- dfDT_codes[,
          (codedTriNames) := ref_window_df(
            .SD[, .SD, .SDcols = just_codes],
            windowSize = window$back,
            windowForward = window$forward,
            binary = binary,
            binaryStanzas = binaryStanzas
          ),
          by = conversations.by,
          .SDcols = initial_cols,
          with = T
      ];
    # }
  }

  ###
  # Convert the generic `V` names to corresponding `adjacency.vector` names
  ###
    colnames(dfDT.co.occurrences)[
      grep("Vd+", colnames(dfDT.co.occurrences))
    ] <- codedTriNames

  ##
  # If units aren't supplied, use all available
  ##
    if (is.null(units.used)) {
      units.used <- dfDT_codes$ENA_UNIT
    }


  ###
  # Trajectory Checks
  ###

  ## Not a Trajectory
  if (enadata$model == "EndPoint") {
    ###
    # Sum each unit found in dfDT.co.occurrences
    ###
    dfDT.summed.units <- dfDT.co.occurrences[ENA_UNIT %in% units.used,lapply(.SD,sum),by=units.by,.SDcols=codedTriNames]
    dfDT.summed.units$ENA_UNIT <- merge_columns_c(dfDT.summed.units, units.by, sep=".");

    enadata$unit.names <- dfDT.summed.units$ENA_UNIT;
  }
  ## Trajectory
  else {
    ## First sum all units within each Trajectory Group (trajectory.by)
    dfDT.summed.traj.by <- dfDT.co.occurrences[
      ENA_UNIT %in% units.used,
      {
        sums <- lapply(.SD, sum)
        data.frame(ENA_ROW_IDX = .GRP, sums); # Return value
      },
      by = c(units.by, trajectory.by),
      .SDcols = (codedTriNames)
    ];
    dfDT.summed.traj.by$ENA_UNIT <- merge_columns_c(
      dfDT.summed.traj.by, units.by, sep = "."
    )
    dfDT.summed.traj.by$TRAJ_UNIT <- merge_columns_c(
      dfDT.summed.traj.by, trajectory.by, sep = "."
    );

    enadata$trajectories$step <- dfDT.summed.traj.by$TRAJ_UNIT;

    # Accumulated
    if (enadata$model == opts$TRAJ_TYPES[1]) {
      dfDT.summed.units <- dfDT.summed.traj.by[
        ENA_UNIT %in% unique(units.used), {
          cols <- colnames(.SD)
          ENA_UNIT <- paste(as.character(.BY), collapse = ".")
          TRAJ_UNIT <- .SD[, c(trajectory.by), with = F]
          inc_cols <- cols[! cols %in% c(trajectory.by, "ENA_ROW_IDX")]
          lag <- ref_window_lag(.SD[, .SD, .SDcols = inc_cols], .N)

          data.table::data.table(
            ENA_ROW_IDX,
            TRAJ_UNIT, lag, ENA_UNIT = ENA_UNIT
          )
        },
        by = c(units.by),
        .SDcols = c(codedTriNames, trajectory.by, "ENA_ROW_IDX")
      ]
      dfDT.summed.units$TRAJ_UNIT <- merge_columns_c(
        dfDT.summed.units, trajectory.by, sep = "."
      )
    }
    # Non-accumulated
    else if (enadata$model == opts$TRAJ_TYPES[2]) {
      dfDT.summed.units <- dfDT.summed.traj.by;
    }
    else {
      stop("Unsupported Model type.");
    }

    dfDT.summed.units$ENA_UNIT <- merge_columns_c(
      dfDT.summed.units, units.by, sep = "."
    )
  }
  ###
  # END: Trajectory Checks
  ###

  ###
  # Name the rows and columns accordingly
  ###
    colnames(dfDT.summed.units)[
      grep("Vd+", colnames(dfDT.summed.units))
    ] <- codedTriNames

  ###
  # Set attributes
  #
  # TODO Most of this should be moved to a more prominent spot on ENAdata
  ###
    adjRows <- triIndices(length(codes)) + 1
    codedRow1 <- codes[adjRows[1, ]]
    codedRow2 <- codes[adjRows[2, ]]
    attr(dfDT.summed.units, "adjacency.matrix") <- rbind(codedRow1, codedRow2)
    attr(dfDT.summed.units, "adjacency.codes") <- codedTriNames
    attr(dfDT.summed.units, opts$UNIT_NAMES)  <- dfDT.summed.units[,
        .SD, with = T, .SDcols = units.by]

    enadata$adjacency.matrix <- rbind(codedRow1, codedRow2)
    enadata$accumulated.adjacency.vectors <- dfDT.co.occurrences
    enadata$adjacency.vectors <- dfDT.summed.units
  ###
  # END: Set attributes
  ###

  return(enadata);
}
"""
lws_positions_sq = """
# Ellipsoidal scaling version
lws.positions.sq <- function(enaset) {
  points = as.matrix(enaset$points)
  weights = as.matrix(enaset$line.weights)
  positions = lws_lsq_positions(weights, points, ncol(points));

  node.positions = positions$nodes;
  rownames(node.positions) = enaset$enadata$codes;

  return(list("node.positions" = node.positions, "centroids" = positions$centroids))
}


lws.positions.sq.R6 <- function(enaset) {
  positions = lws_lsq_positions(enaset$line.weights, enaset$points.rotated, enaset$get("dimensions"));

  node.positions = positions$nodes;
  rownames(node.positions) = enaset$enadata$codes;

  return(list("node.positions" = node.positions, "centroids" = positions$centroids))

}
"""
cohens_d = """
##
#' Cohen's d calculation
#'
#' @title Cohen's d
#'
#' @description Calculate Conhen's d
#'
#' @details [TBD]
#'
#' @param x [TBD]
#' @param y [TBD]
#'
#' @export
#' @return numeric Cohen's d calculation
fun_cohens.d <- function(x, y) {
  lx <- length(x)- 1
  ly <- length(y)- 1
  md  <- abs(mean(x) - mean(y))        ## mean difference (numerator)
  csd <- lx * var(x) + ly * var(y)
  csd <- csd/(lx + ly)
  csd <- sqrt(csd)                     ## common sd computation

  cd  <- md/csd
  return(cd)## cohen's d
}
"""
connection_matrix = """
#' Connection counts as square matrix
#'
#' @param x ena.set or ena.connections (i.e. set$connection.counts)
#'
#' @return matrix
#' @export
connection.matrix <- function(x) {
  if(is(x, "ena.set")) {
    connections <- x$connection.counts
  } else {
    connections <- x
  }
  if(!is(connections, "ena.connections")) {
    stop("Unable to find connections. `x` must be connections from an ena.set or an ena.set")
  }

  simplify <-  (nrow(connections) == 1)
  cm <- as.matrix(connections, square = T, simplify = simplify)
  if(simplify == FALSE && is.list(cm))
    names(cm) <- connections$ENA_UNIT

  return(cm);
}
"""
ena = """
#####
#' @title Wrapper to generate, and optionally plot, an ENA model
#'
#' @description Generates an ENA model by constructing a dimensional reduction
#' of adjacency (co-occurrence) vectors as defined by the supplied
#' conversations, units, and codes.
#'
#' @details This function generates an ena.set object given a data.frame, units,
#' conversations, and codes. After accumulating the adjacency (co-occurrence)
#' vectors, computes a dimensional reduction (projection), and calculates node
#' positions in the projected ENA space. Returns location of the units in the
#' projected space, as well as locations for node positions, and normalized
#' adjacency (co-occurrence) vectors to construct network graphs. Includes options
#' for returning statistical tests between groups of units, as well as plots of units,
#' groups, and networks.
#'
#'
#'
#'
#'
#' @param data data.frame with containing metadata and coded columns
#' @param codes vector, numeric or character, of columns with codes
#' @param units vector, numeric or character, of columns representing units
#' @param conversation  vector, numeric or character, of columns to segment conversations by
#' @param metadata  vector, numeric or character, of columns with additional meta information for units
#' @param model character: EndPoint (default), AccumulatedTrajectory, SeparateTrajectory
#' @param weight.by "binary" is default, can supply a function to call (e.g. sum)
#' @param window MovingStanzaWindow (default) or Conversation
#' @param window.size.back Number of lines in the stanza window (default: 1)
#' @param include.meta [TBD]
#' @param groupVar vector, character, of column name containing group identifiers.
#' If column contains at least two unique values, will generate model using a means rotation (a dimensional reduction maximizing the variance between the means of the two groups)
#' @param groups vector, character, of values of groupVar column used for means rotation, plotting, or statistical tests
#' @param runTest logical, TRUE will run a Student's t-Test and a Wilcoxon test for groups defined by the groups argument
#' @param points logical, TRUE will plot points (default: FALSE)
#' @param mean logical, TRUE will plot the mean position of the groups defined in the groups argument (default: FALSE)
#' @param network logical, TRUE will plot networks (default: TRUE)
#' @param networkMultiplier numeric, scaling factor for non-subtracted networks (default: 1)
#' @param subtractionMultiplier numeric, scaling factor for subtracted networks (default: 1)
#' @param unit vector, character, name of a single unit to plot
#' @param include.plots logical, TRUE will generate plots based on the model (default: TRUE)
#' @param print.plots logical, TRUE will show plots in the Viewer(default: FALSE)
#' @param ... Additional parameters passed to set creation and plotting functions
#'
#' @examples
#' data(RS.data)
#'
#' rs = ena(
#'   data = RS.data,
#'   units = c("UserName","Condition", "GroupName"),
#'   conversation = c("Condition","GroupName"),
#'   codes = c('Data',
#'             'Technical.Constraints',
#'             'Performance.Parameters',
#'             'Client.and.Consultant.Requests',
#'             'Design.Reasoning',
#'             'Collaboration'),
#'   window.size.back = 4,
#'   print.plots = FALSE,
#'   groupVar = "Condition",
#'   groups = c("FirstGame", "SecondGame")
#' )
#'
#' @return ena.set object
#' @export
#####
ena <- function(
  data,
  codes,
  units,
  conversation,
  metadata = NULL,
  model = c("EndPoint", "AccumulatedTrajectory", "SeparateTrajectory"),
  weight.by = "binary",
  window = c("MovingStanzaWindow", "Conversation"),
  window.size.back = 1,
  include.meta = TRUE,
  groupVar = NULL,
  groups = NULL,
  runTest = FALSE,
  points = FALSE,
  mean = FALSE,
  network = TRUE,
  networkMultiplier = 1,
  subtractionMultiplier = 1,
  unit = NULL,
  include.plots = T,
  print.plots = F,
  ...
) {
  set <- ena.set.creator(
    data = data,
    codes = codes,
    units = units,
    conversation = conversation,
    metadata = metadata,
    model = model,
    weight.by = weight.by,
    window = window,
    window.size.back = window.size.back,
    include.meta = include.meta,
    groupVar = groupVar,
    groups = groups,
    runTest = runTest,
    # testType = testType,
    ...
  )

  if (include.plots) {
    set <- ena.plotter(
      set = set,
      groupVar = groupVar,
      groups = groups,
      points = points,
      mean = mean,
      network = network,
      networkMultiplier = networkMultiplier,
      subtractionMultiplier = subtractionMultiplier,
      unit = unit,
      print.plots = print.plots,
      ...
    )
  }

  return(set)
}
"""
ena_accumulate_data = """
##
#' @title Accumulate data from a data frame into a set of adjacency (co-occurrence) vectors
#'
#' @description This function initializes an ENtimea object, processing conversations from coded data to generate adjacency (co-occurrence) vectors
#'
#' @details ENAData objects are created using this function. This accumulation receives
#' separate data frames for units, codes, conversation, and optionally, metadata. It
#' iterates through the data to create an adjacency (co-occurrence) vector corresponding
#' to each unit - or in a trajectory model multiple adjacency (co-occurrence) vectors for
#' each unit.
#'
#' In the default MovingStanzaWindow model, co-occurrences between codes are
#' calculated for each line k in the data between line k and the window.size.back-1 previous
#' lines and window.size.forward-1 subsequent lines in the same conversation as line k.
#'
#' In the Conversation model, co-occurrences between codes are calculated across all lines in
#' each conversation. Adjacency (co-occurrence) vectors are constructed for each unit u by
#' summing the co-occurrences for the lines that correspond to u.
#'
#' Options for how the data is accumulated are endpoint, which produces one adjacency (co-occurrence)
#' vector for each until summing the co-occurrences for all lines, and two trajectory models:
#' AccumulatedTrajectory and SeparateTrajectory. Trajectory models produce an adjacency
#' (co-occurrence) model for each conversation for each unit. In a SeparateTrajectory model,
#' each conversation is modeled as a separate network. In an AccumulatedTrajectory model, the
#' adjacency (co-occurrence) vector for the current conversation includes the co-occurrences
#' from all previous conversations in the data.
#'
#' @export
#'
#' @param units A data frame where the columns are the properties by which units will be identified
#' @param conversation A data frame where the columns are the properties by which conversations will be identified
#' @param codes A data frame where the columns are the codes used to create adjacency (co-occurrence) vectors
#' @param metadata (optional) A data frame with additional columns of metadata to be associated with each unit in the data
#' @param model A character, choices: EndPoint (or E), AccumulatedTrajectory (or A), or SeparateTrajectory (or S); default: EndPoint. Determines the ENA model to be constructed
#' @param weight.by (optional) A function to apply to values after accumulation
#' @param mask (optional) A binary matrix of size ncol(codes) x ncol(codes). 0s in the mask matrix row i column j indicates that co-occurrence will not be modeled between code i and code j
#' @param window A character, choices are Conversation (or C), MovingStanzaWindow (MSW, MS); default MovingStanzaWindow. Determines how stanzas are constructed, which defines how co-occurrences are modeled
#' @param window.size.back A positive integer, Inf, or character (INF or Infinite), default: 1. Determines, for each line in the data frame, the number of previous lines in a conversation to include in the stanza window, which defines how co-occurrences are modeled
#' @param window.size.forward (optional) A positive integer, Inf, or character (INF or Infinite), default: 0. Determines, for each line in the data frame, the number of subsequent lines in a conversation to include in the stanza window, which defines how co-occurrences are modeled
#' @param ... additional parameters addressed in inner function
#' @param include.meta Locigal indicating if unit metadata should be attached to the resulting ENAdata object, default is TRUE
#' @param as.list R6 objects will be deprecated, but if this is TRUE, the original R6 object will be returned, otherwise a list with class `ena.set`
#'
#' @seealso \code{\link{ENAdata}}, \code{\link{ena.make.set}}
#'
#' @return \code{\link{ENAdata}} object with data [adjacency (co-occurrence) vectors] accumulated from the provided data frames.
#'
##
ena.accumulate.data <- function(
  units = NULL,
  conversation = NULL,
  codes = NULL,
  metadata = NULL,
  model = c("EndPoint", "AccumulatedTrajectory", "SeparateTrajectory"),
  weight.by = "binary",
  window = c("MovingStanzaWindow", "Conversation"),
  window.size.back = 1,
  window.size.forward = 0,
  mask = NULL,
  include.meta = T,
  as.list = T,
  ...
) {
  if (is.null(units) || is.null(conversation) || is.null(codes)) {
    stop("Accumulation requires: units, conversation, and codes");
  }
  if (nrow(units) != nrow(conversation) || nrow(conversation) != nrow(codes)) {
    stop("Data Frames do not have the same number of rows");
  }

  df <- cbind(units, conversation);
  df <- cbind(df, codes);

  metadata <- data.table::as.data.table(metadata)
  if (!is.null(metadata) && nrow(metadata) == nrow(df)) {
    df <- cbind(df, metadata);
  }

  model <- match.arg(model)
  window <- match.arg(window)

  units.by <- colnames(units);
  conversations.by <- colnames(conversation);
  if (identical(window, "Conversation")) {
    conversations.by <- c(conversations.by, units.by);
    window.size.back <- window;
  }
  else if (identical(window, "MovingStanzaWindow")) {
    if( grepl(pattern = "inf", x = window.size.back, ignore.case = T)) {
      window.size.back <- Inf
    }
    if( grepl(pattern = "inf", x = window.size.forward, ignore.case = T)) {
      window.size.forward <- Inf
    }
  }

  data <- ENAdata$new(
    file = df,
    units = units,
    units.by = units.by,
    conversations.by = conversations.by,
    codes = codes,
    window.size.back = window.size.back,
    window.size.forward = window.size.forward,
    weight.by = weight.by,
    model = model,
    mask = mask,
    include.meta = include.meta,
    ...
  );
  data$process()

  data$function.call <- sys.call()

  if(as.list) {
    data <- ena.set(data)
  } else {
    warning(paste0("Usage of R6 data objects is deprecated and may be removed ",
      "entirely in a future version. Consider upgrading to the new data ",
      "object."))
  }

  data
}
"""
ena_accumulate_data_file = """
##
# @title Accumulate Data from csv
#
# @description This function accumulates rows of data.
#
# @details [TBD]
#
#@export
#
# @param file The csv file location or data.frame for the function
# @param units.used Delimits columns based on the units (which specific units to use)
# @param units.by unit columns to accumulate by
# @param conversations.by Columns used in the conversation
# @param codes Columns used based on codes
# @param window.size.back Number of lines back to include window in stanza
# @param window.size.forward Number of lines forward in stanza window
# @param binary [TBD]
# @param model [TBD]
# @param window [TBD]
# @param weight.by [TBD]
# @param binary.stanzas [TBD]
# @param mask [TBD]
# @param ... additional parameters addressed in inner function
#
#
# @seealso \code{\link{ena.make.set}}
#
# @examples
# \dontrun{
# codeNames = c(
#   "E.data","S.data","E.design","S.design","S.professional","E.client",
#   "V.client","E.consultant","V.consultant","S.collaboration","I.engineer",
#   "I.intern","K.actuator","K.rom","K.materials","K.power"
# )
#
# df.file <- system.file("extdata", "rs.data.csv", package="rENA")
#
# # Given a csv file location
# ena.accumulate.data(
#   df.file, units.by = c("UserName","Condition"),
#   conversations.by = c("ActivityNumber","GroupName"),
#   codes = codeNames
# )
# }
# @return \code{\link{ENAdata}} class object with accumulated data
#
##
ena.accumulate.data.file <- function(
  file,
  units.used = NULL,
  conversations.used = NULL,
  units.by,
  conversations.by,
  codes = NULL,
  model = c("EndPoint",
            "AccumulatedTrajectory",
            "SeparateTrajectory"),
  window = c("Moving Stanza", "Conversation"),
  window.size.back = 1,
  window.size.forward = 0,
  weight.by = "binary",
  binary.stanzas = F,
  mask = NULL,
  include.meta = T,
  as.list = T,
  ...
) {
  if(is.null(file) ||
     is.null(units.by) ||
     is.null(conversations.by) || is.null(codes)
  ) {
    stop("Accumulation: file, units.by, conversations.by, and codes")
  }

  units <- NULL;
  model <- match.arg(model);
  window <- match.arg(window);

  if (identical(window, "Conversation")) {
    conversations.by = c(conversations.by, units.by);
    window.size.back = window;
  }
  data = ENAdata$new(
    file = file,
    units = units,
    units.used = units.used,
    units.by = units.by,
    conversations.by = conversations.by,
    codes = codes,
    window.size.back = window.size.back,
    window.size.forward = window.size.forward,
    weight.by = weight.by,
    model = model,
    mask = mask,
    include.meta = include.meta,
    ...
  );
  data$process();

  data$function.call = sys.call();
  # output = match.arg(output);
  # if(output == "json") {
  #   output.class = get(class(data))
  #
  #   if(is.null(output.fields)) {
  #     output.fields = names(output.class$public_fields)
  #   }
  #
  #   r6.to.json(data, o.class = output.class, o.fields = output.fields)
  # }
  #else

  if(as.list) {
    data = ena.set(data);
  } else {
    warning("Usage of R6 data objects is deprecated and may be removed entirely in a future version. Consider upgrading to the new data object.")
  }
  data
}
"""
ena_conversations = """
##
#' @title Find conversations by unit
#'
#' @description Find rows of conversations by unit
#'
#' @details [TBD]
#'
#' @param set [TBD]
#' @param units [TBD]
#' @param units.by [TBD]
#' @param codes [TBD]
#' @param conversation.by [TBD]
#' @param window [TBD]
#' @param conversation.exclude [TBD]
#'
#' @examples
#' data(RS.data)
#'
#' codeNames = c('Data','Technical.Constraints','Performance.Parameters',
#'               'Client.and.Consultant.Requests','Design.Reasoning',
#'               'Collaboration');
#'
#' accum = ena.accumulate.data(
#'   units = RS.data[,c("Condition","UserName")],
#'   conversation = RS.data[,c("Condition","GroupName")],
#'   metadata = RS.data[,c("CONFIDENCE.Change","CONFIDENCE.Pre",
#'                         "CONFIDENCE.Post","C.Change")],
#'   codes = RS.data[,codeNames],
#'   model = "EndPoint",
#'   window.size.back = 4
#' );
#' set = ena.make.set(
#'   enadata = accum,
#'   rotation.by = ena.rotate.by.mean,
#'   rotation.params = list(accum$meta.data$Condition=="FirstGame",
#'                          accum$meta.data$Condition=="SecondGame")
#' );
#' ena.conversations(set = RS.data,
#'   units = c("FirstGame.steven z"), units.by=c("Condition","UserName"),
#'   conversation.by = c("Condition","GroupName"),
#'   codes=codeNames, window = 4
#' )
#'
#' @export
#' @return list containing row indices representing conversations
##
ena.conversations = function(set, units, units.by=NULL, codes=NULL, conversation.by = NULL, window = 4, conversation.exclude = c()) {
  # rawData = data.table::copy(set$enadata$raw);
  if(is.null(units.by)) {
    units.by = set$`_function.params`$units.by;
  }
  # conversation.by = set$enadata$function.params$conversations.by;
  # window = set$enadata$function.params$window.size.back;
  # rawAcc = data.table::copy(set$enadata$accumulated.adjacency.vectors);
  if(is(set, "ena.set")) {
    rawAcc2 = set$model$raw.input
  } else {
    rawAcc2 = data.table::data.table(set) #$enadata$raw);
  }

  # rawAcc$KEYCOL = merge_columns_c(rawAcc, conversation.by)
  rawAcc2$KEYCOL = merge_columns_c(rawAcc2, conversation.by)

  # conversationsTable = rawAcc[, paste(.I, collapse = ","), by = c(conversation.by)]
  conversationsTable2 = rawAcc2[, paste(.I, collapse = ","), by = c(conversation.by)]

  # rows = sapply(conversationsTable$V1, function(x) as.numeric(unlist(strsplit(x, split=","))),USE.NAMES = T)
  rows2 = lapply(conversationsTable2$V1, function(x) as.numeric(unlist(strsplit(x, split=","))))
  # browser()
  # names(rows) = merge_columns_c(conversationsTable,conversation.by); #unique(rawAcc[,KEYCOL])
  names(rows2) = merge_columns_c(conversationsTable2,conversation.by); #unique(rawAcc[,KEYCOL])

  # unitRows = merge_columns_c(rawAcc[,c(units.by),with=F], units.by)
  unitRows2 = merge_columns_c(rawAcc2[,c(units.by),with=F], units.by)

  # adjCol = set$enadata$adjacency.matrix[1,] %in%  codes[1] & set$enadata$adjacency.matrix[2,] %in% codes[2]
  # adjColName = paste("adjacency.code.", which(adjCol), sep = "")
  # codedUnitRows = which(unitRows %in% units & rawAcc[[adjColName]] == 1)

  codedRows = rawAcc2[, rowSums(.SD), .SDcols = codes] > 0
  codedUnitRows2 = which(unitRows2 %in% units & codedRows)
  codedUnitRows2 = codedUnitRows2[!(codedUnitRows2 %in% as.vector(unlist(rows2[conversation.exclude])))]
  # codedUnitRowConvs = rawAcc[codedUnitRows,KEYCOL];
  codedUnitRowConvs2 = rawAcc2[codedUnitRows2,KEYCOL];

  # codedUnitRowConvsAll = NULL;
  codedUnitRowConvsAll2 = NULL;
  unitRowsNotCooccurred = c()
  if(length(codedUnitRows2) > 0) {
    codedUnitRowConvsAll = unique(unlist(sapply(X = 1:length(codedUnitRows2), simplify = F, FUN = function(x) {
      thisConvRows = rows2[[codedUnitRowConvs2[x]]]
      thisRowInConv = which(thisConvRows == codedUnitRows2[x])
      winUse = ifelse(is.infinite(window), thisRowInConv, window)
      thisRowAndWindow = rep(thisRowInConv,winUse) - (winUse-1):0
      coOccursFound = all(rawAcc2[thisConvRows[thisRowAndWindow[thisRowAndWindow > 0]], lapply(.SD, sum), .SDcols=codes] > 0)
      if(coOccursFound) {
        thisConvRows[thisRowAndWindow[thisRowAndWindow > 0]]
      } else {
        unitRowsNotCooccurred <<- c(unitRowsNotCooccurred, thisConvRows[thisRowInConv])
        # coOccursFound
        NULL
      }
    })))
  }
  return(list(
    conversations = as.list(rows2),
    unitConvs = unique(rawAcc2[codedUnitRows2,KEYCOL]),
    allRows = codedUnitRowConvsAll,
    unitRows = codedUnitRows2,
    toRemove = unitRowsNotCooccurred
  ));
}

"""
ena_correlations = """
###
#' Calculate the correlations
#'
#' @description Calculate both Spearman and Pearson correlations for the
#' provided ENAset
#'
#' @param enaset ENAset to run correlations on
#' @param dims The dimensions to calculate the correlations for. Default: c(1,2)
#'
#' @return Matrix of 2 columns, one for each correlation method, with the corresponding
#' correlations per dimension as the rows.
#'
#' @export
###
ena.correlations <- function(enaset, dims = c(1:2)) {
  pComb = combn(nrow(enaset$points),2)
  point1 = pComb[1,]
  point2 = pComb[2,]

  points = as.matrix(enaset$points)
  centroids = as.matrix(enaset$model$centroids)
  svdDiff = matrix(points[point1, dims] - points[point2, dims], ncol=length(dims), nrow=length(point1))
  optDiff = matrix(centroids[point1, dims] - centroids[point2, dims], ncol=length(dims), nrow=length(point1))

  correlations = as.data.frame(mapply(function(method) {
    lapply(dims, function(dim) {
      cor(as.numeric(svdDiff[,dim]), as.numeric(optDiff[,dim]), method=method)
    });
  }, c("pearson","spearman")))

  return(correlations);
}
"""
ena_generate = """
# # ##
# # # @title Accumulate and Generate
# # #
# # # @description Accumulate and Generate
# # #
# # # @details [TBD]
# # #
# # # @param file [TBD]
# # # @param window.size.back [TBD]
# # # @param units.by [TBD]
# # # @param conversations.by [TBD]
# # # @param code [TBD]
# # # @param units.used [TBD]
# # # @export
# # # @return list containing the accumulation and set
# # ##
# # ena.generate <- function(
# #   file,
# #   window.size.back,
# #   units.by,
# #   conversations.by,
# #   code,
# #   scale.nodes = T,
# #   units.used = NULL,
# #   ...
# # ) {
# #   args = list(...);
# #   conversations.used = NULL;
# #   weight.by = "binary";
# #   if(!is.null(args$conversations.used)) {
# #     conversations.used = args$conversations.used
# #     file$KEYCOL = merge_columns_c(file,conversations.by)
# #     file = file[file$KEYCOL %in% conversations.used,]
# #   }
# #   if(!is.null(args$weight.by)) {
# #     weight.by = args$weight.by
# #   }
# #
# #   accum = ena.accumulate.data.file(
# #     file = file,
# #     window.size.back = window.size.back,
# #     units.by = make.names(units.by),
# #     units.used = units.used,
# #     model = "EndPoint",
# #     conversations.by = make.names(conversations.by),
# #     codes = make.names(code),
# #     ...
# #   )
# #
# #   rotate.groups = NULL
# #   if(!is.null(args$rotate.by)) {
# #     rotate.meta = accum$metadata[accum$metadata$ENA_UNIT %in% accum$unit.names,]
# #     rotate.col = accum$metadata[accum$metadata$ENA_UNIT %in% accum$unit.names,][[names(args$rotate.by)[1]]]
# #     rotate.groups = list(
# #       rotate.col == args$rotate.by[[1]][1],
# #       rotate.col == args$rotate.by[[1]][2]
# #     )
# #   }
# #   set = ena.make.set(
# #     enadata = accum,
# #     norm.by = ifelse((is.null(args$sphere.norm) || args$sphere.norm==T),sphere_norm_c,dont_sphere_norm_c),
# #     rotation.by = if(is.null(rotate.groups)) ena.svd else ena.rotate.by.mean, #ifelse(is.null(rotate.groups), NULL, ena.rotate.by.mean),
# #     rotation.params = rotate.groups,
# #     ...
# #   )
# #
# # @param file [TBD]
# # @param window.size.back [TBD]
# # @param units.by [TBD]
# # @param conversations.by [TBD]
# # @param code [TBD]
# # @param units.used [TBD]
# # @export
# # @return list containing the accumulation and set
# ##
# ena.generate <- function(
#   file,
#   window.size.back,
#   units.by,
#   conversations.by,
#   code,
#   scale.nodes = T,
#   units.used = NULL,
#   dimensions = 6,
#   include.meta = F,
#   ...
# ) {
#   startedTime = as.numeric(Sys.time())
#   args = list(...);
#   unit.groups = NULL;
#   conversations.used = NULL;
#   weight.by = "binary";

#   if(!is.null(args$conversations.used)) {
#     conversations.used = args$conversations.used
#     file$KEYCOL = rENA:::merge_columns_c(file, make.names(conversations.by))
#     file = file[file$KEYCOL %in% conversations.used,]
#   }

#   if(!is.null(args$weight.by)) {
#     weight.by = args$weight.by;
#   }
#   if(!is.null(args$unit.groups)){
#     if(is.data.frame(args$unit.groups)) {
#       unit.groups = args$unit.groups$units;
#       names(unit.groups) = args$unit.groups$name;
#     } else {
#       unit.groups = list();
#       group.json = jsonlite::fromJSON(args$unit.groups)
#       for(grp in 1:length(group.json$name)) {
#         unit.groups[group.json$name[grp]] = group.json$units[grp];
#       }
#     }
#   }

#   accum = ena.accumulate.data.file(
#     file = file,
#     window.size.back = window.size.back,
#     units.by = make.names(units.by),
#     units.used = units.used,
#     model = "EndPoint",
#     conversations.by = make.names(conversations.by),
#     codes = make.names(code),
#     include.meta = include.meta,
#     ...
#   )

#   rotate.groups = NULL
#   if(!is.null(args$rotate.by)) {
#     # rotate.meta = accum$metadata[accum$metadata$ENA_UNIT %in% accum$unit.names,]
#     # rotate.col = accum$metadata[accum$metadata$ENA_UNIT %in% accum$unit.names,][[make.names(names(args$rotate.by)[1])]]
#     # rotate.groups = list(
#     #   rotate.col == args$rotate.by[[1]][1],
#     #   rotate.col == args$rotate.by[[1]][2]
#     # )
#     rotate.groups = lapply(args$rotate.by, function(x) accum$unit.names %in% x )
#   }

#   use.to.norm = rENA:::fun_skip_sphere_norm;
#   if(is.null(args$sphere.norm) || args$sphere.norm == T) {
#     use.to.norm = rENA:::fun_sphere_norm
#   }
#   rotation.set = NULL
#   if(!is.null(args$rotation.matrix)) {
#     rotation.set = ENARotationSet$new(
#       rotation = args$rotation.matrix$rotation$rotation,
#       node.positions = args$rotation.matrix$rotation$node.positions,
#       codes = args$rotation.matrix$codes
#     );
#     colnames(rotation.set$rotation) = args$rotation.matrix$dimensions;
#   }
#   set = ena.make.set(
#     enadata = accum,
#     norm.by = use.to.norm,
#     rotation.by = if(is.null(rotate.groups)) rENA:::ena.svd else rENA:::ena.rotate.by.mean,
#     rotation.params = rotate.groups,
#     rotation.set = rotation.set,
#     dimensions = dimensions,
#     ...
#   )

#   use.dimensions = 1:2
#   # if(!is.null(args$keep.dimensions)) {
#   #   use.dimensions = which(colnames(set$points.rotated) %in% args$keep.dimensions)
#   # }

#   group.names = NULL;
#   if(length(units.by)>1) {
#     group.names = unique(set$enadata$units[[make.names(units.by)[[1]]]])
#   } else {
#     group.names = units.by
#   }
#   group.cnt = length(group.names);
#   conf.ints = list();
#   outlier.ints = matrix(0, nrow=(group.cnt), ncol=(2));

#   set$points.rotated.scaled = set$points.rotated;
#   scaleFactor = 1.0
#   if(scale.nodes == T) {
#     np.min.x = min(set$node.positions[,use.dimensions[1]])
#     np.min.y = min(set$node.positions[,use.dimensions[2]])
#     rp.min.x = min(set$points.rotated[,use.dimensions[1]])
#     rp.min.y = min(set$points.rotated[,use.dimensions[2]])
#     maxMin = abs(max(np.min.x / rp.min.x, np.min.y / rp.min.y))

#     np.max.x = max(set$node.positions[,use.dimensions[1]])
#     np.max.y = max(set$node.positions[,use.dimensions[2]])
#     rp.max.x = max(set$points.rotated[,use.dimensions[1]])
#     rp.max.y = max(set$points.rotated[,use.dimensions[2]])
#     maxMax = abs(max(np.max.x / rp.max.x, np.max.y / rp.max.y))
#     scaleFactor = min(maxMin, maxMax)
#     # set$points.rotated = set$points.rotated * scaleFactor;
#     set$points.rotated.scaled = set$points.rotated * scaleFactor;
#   }

#   # groups = NULL
#   group.method = "mean"
#   if(!is.null(args$weight.network.by) && (args$weight.network.by %in% c("mean","sum"))) {
#     group.method = args$weight.network.by;
#   }

#   groups = list()
#   if(!is.null(unit.groups) && length(unit.groups) > 0){
#     # for(i in 1:length(names(unit.groups))) {
#     #   groups[[length(groups)+1]] = ena.unit.group(set, set$enadata$unit.names[set$enadata$unit.names %in% unit.groups[[i]]], name = names(unit.groups)[i], method = group.method, scaleFactor = scaleFactor)
#     # }
#     groups = lapply(names(unit.groups), function(nm) {
#       ena.unit.group(
#         set,
#         set$enadata$unit.names[set$enadata$unit.names %in% unit.groups[[nm]]],
#         name = nm,
#         method = group.method,
#         scaleFactor = scaleFactor
#         # ,keep.dimensions = use.dimensions
#       )
#     })
#   }

#   if(
#     !is.null(args$output) && args$output == "save" &&
#     !is.null(args$output.to)
#   ) {
#     setName = tools::file_path_sans_ext(basename(args$output.to))
#     env = environment()
#     assign(x = setName, value = set, envir = env);
#     env[[setName]] = get(x = setName, envir = env)

#     tmp <- tempfile(fileext = ".rdata")
#     on.exit(unlink(tmp))
#     save(list = c(setName), file = tmp, envir = env)
#     bucket <- aws.s3::get_bucketname(args$output.to)
#     object <- aws.s3:::get_objectkey.character(args$output.to)
#     return(aws.s3::put_object(file = tmp, bucket = bucket, object = object));
#   }
#   else {
#     nodes = data.frame(set$node.positions);
#     nodes$weight = rep(0, nrow(nodes))
#     node.rows = rownames(set$node.positions);
#     estimate.over.units = (!(set$enadata$unit.names %in% args$units.exclude))
#     weights = matrix(0, ncol=nrow(set$node.positions), nrow=length(which(estimate.over.units)));

#     colnames(weights) = node.rows
#     network.scaled = set$line.weights[estimate.over.units,];
#     # if(!is.null(scale.weights) && scale.weights == T) {
#     #   network.scaled = network.scaled * (1 / max(abs(network.scaled)));
#     # }

#     mat = set$enadata$adjacency.matrix;
#     # for (x in 1:nrow(network.scaled)) {
#     #   weights[x, ] = sapply(node.rows, function(y) {
#     #     # sum(network.scaled[x,as.logical(colSums(!is.na(apply(mat,2,match, y))))])
#     #     sum(network.scaled[x, as.logical(colSums(mat == y))])
#     #   })
#     #   # network.thickness = network.scaled[x,] #scales::rescale(abs(network.scaled[x,]), thickness);
#     #   # for (i in 1:ncol(mat)) {
#     #   #   weights[x,node.rows==mat[1,i]] = weights[x,node.rows==mat[1,i]] + network.thickness[i];
#     #   # }
#     # }
#     weights = sapply(node.rows, function(x) rowSums(network.scaled[,as.logical(colSums(mat == x) )]))

#     # #weights = t(apply(weights, 1, scales::rescale, c(1,ncol(weights))));
#     weights = scales::rescale(weights, c(1,ncol(weights)));

#     set$line.weights[estimate.over.units,] = set$line.weights[estimate.over.units,];

#     # If not included, remove the weights as to not effect the scaling
#     set$line.weights[!estimate.over.units,] = 0
#     scaleRange = c(min(set$line.weights[estimate.over.units,]) ,1);
#     if(scaleRange[1] < 0.1 && min(set$line.weights)>0) {
#       scaleRange[1] = 0.1;
#     }

#     set$line.weights = scales::rescale(set$line.weights, to=scaleRange, from=range(set$line.weights, na.rm = T, finite = T))
#     # adjRows = triIndices(length(code)) + 1
#     # codedRow1 = code[adjRows[1,]];
#     # codedRow2 = code[adjRows[2,]];

#     tmp = getwd();
#     sess = regexec("temp/(x[^/]*)/workspace", tmp)[[1]]
#     assign("set", set, envir = parent.frame())

#     dimension.names = paste("SVD",1:ncol(set$points.rotated), sep="")
#     if(length(set$function.params$rotation.params) == 2) dimension.names[1] = "MR1"
#     # if(length(args$plotted.nodes) == 2) {
#     #   methods = ena.methods(enaset = set, tool = "webENA", tool.version = "0.1.0", comparison = "parametric", comparison.groups = args$plotted.nodes)
#     # } else {
#     #   methods = ena.methods(enaset = set, tool = "webENA", tool.version = "0.1.0")
#     # }
#     doneTime = as.numeric(Sys.time())

#     ### Limit dimensions
#     # set$points.rotated = set$points.rotated[,use.dimensions]
#     # set$points.rotated.scaled = set$points.rotated.scaled[,use.dimensions]
#     # set$node.positions = set$node.positions[,use.dimensions]
#     # set$rotation.set$rotation = set$rotation.set$rotation[,use.dimensions]
#     return(list(
#       codes = make.names(code),
#       adjacency.matrix = mat, #rbind(codedRow1, codedRow2),
#       set = set,
#       # methods = readChar(methods, file.info(methods)$size),
#       custom.rotation = if(!is.null(rotation.set)) T else F,
#       custom.rotation.set = rotation.set,
#       dimensions = dimension.names, #colnames(set$points.rotated),
#       session = substr(tmp, start=sess[2], stop=sess[2]+attr(sess, "match.length")[2]-1),
#       # groups = groups,
#       groups = groups,
#       scaled = scale.nodes,
#       node.sizes = weights,
#       esitmated.over = args$units.exclude,
#       edge.saturation = scales::rescale(set$line.weights, c(0.25,1)),
#       edge.opacity = scales::rescale(set$line.weights, c(0.3,1)),
#       startedTime = startedTime,
#       doneTime = doneTime,
#       durationTime = doneTime - startedTime
#     ));
#   }
# }
"""
ena_get_stats = """
##
# @title Generate ENA Set
# @description Generate an ENA set from a givent ENA data object
#
#
#
# @export
##
ena.get.stats <- function(

) {

}
"""
ena_group = """
##
#' @title Compute summary statistic for groupings of units using given method (typically, mean)
#'
#' @description Computes summary statistics for groupings (given as vector) of units in ena data using given method (typically, mean); computes summary statistic for point locations and edge weights for each grouping
#'
#' @export
#'
#' @param enaset An \code{\link{ENAset}} or a vector of values to group.
#' @param by A vector of values the same length as units. Uses rotated points for group positions and normed data to get the group edge weights
#' @param method A function that is used on grouped points. Default: mean().  If `enaset` is an ENAset, enaset$points.rotated will be groups using `mean` regardless of `method` provided
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
#' means = ena.group(set, by=accum$metadata$Condition)
#'
#'
#' @return A list containing names, points, and edge weights for each of the unique groups formed by the function
##
ena.group <- function(
  enaset = NULL,   #ENAset object to form groups from
  by = NULL, #Vector of values  the same length as units.
  method = mean  #method by which to form groups from specified attribute/vector of values
) {
  run.method = function(pts, m = method) {
    points.dt = data.table::data.table(pts);
    if(is.logical(by)) {
      points.dt.means = points.dt[by, lapply(.SD, m),]; # by=by];
    } else if(all(by %in% colnames(pts))) {
      points.dt.means = points.dt[, lapply(.SD, m), by=by];
    } else {
      points.dt.means = as.data.frame(aggregate(points.dt, by = list(by), FUN = m)) #"mean"))
      rownames(points.dt.means) = points.dt.means$Group.1
      points.dt.means = points.dt.means[,colnames(points.dt)]
      # agg.df[as.vector(unique(group.by)),]u
      # return (points.dt.means[as.vector(unique(by)),]);
      points.dt.means[['ENA_GROUP_NAME']] = rownames(points.dt.means)
      return(points.dt.means[which(rownames(points.dt.means) %in% unique(by)),])
    }
    return(as.data.frame(points.dt.means[,colnames(points.dt),with=F]))
  }

  if(is.character(method)) {
    method = get(method)
  }

  if("ENAset" %in% class(enaset)) {
    return(list(
      "names" = as.vector(unique(by)),
      "points" = run.method(enaset$points.rotated, m = mean),
      "line.weights" = run.method(enaset$line.weights)
    ));
  } else {
    return(run.method(enaset))
  }
}
"""
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
            "expected objects:",
            "Expected: c('node.positions','centroids')",
            "Received: ", names(positions), sep = ""))
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
ena_optimize_set = """
##
# @title Generate ENA Set
# @description Generate an ENA set from a givent ENA data object
#
#
#
# @export
##
ena.optimize.set <- function(

) {

}
"""
ena_plot = """
##
#' @title Generate a plot of an ENAset
#'
#' @description Generates an a plot from a given ENA set object
#'
#' @details This function defines the axes and other features of a plot for displaying an ENAset; generates an ENAplot object that can used to plot points, network graphs, and other information from an ENAset
#'
#' @export
#'
#' @param enaset The \code{\link{ENAset}} that will be used to generate a plot
#' @param title A character used for the title of the plot, default: ENA Plot
#' @param dimension.labels A character vector containing labels for the axes, default: c(X, Y)
#' @param font.size An integer determining the font size for graph labels, default: 10
#' @param font.color A character determining the color of label font, default: black
#' @param font.family A character determining the font type, choices: Arial, Courier New, Times New Roman, default: Arial
#' @param scale.to "network" (default), "points", or a list with x and y ranges. Network and points both scale to the c(-max, max) of the corresponding data.frame
#' @param ... additional parameters addressed in inner function
#'
#'
#' @seealso \code{\link{ena.make.set}}, \code{\link{ena.plot.points}}
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
#' plot = ena.plot(set)
#'
#' group1.points = set$points.rotated[set$enadata$units$Condition == "FirstGame",]
#' plot = ena.plot.points(plot, points = group1.points);
#' print(plot);
#'
#' @return \code{\link{ENAplot}} used for plotting an ENAset

##
ena.plot <- function(
  enaset,

  title = "ENA Plot",

  dimension.labels = c("",""),

  font.size = 10,
  font.color = "#000000",
  font.family = c("Arial", "Courier New"),
  scale.to = "network", #, "points"),
  ...
) {
  if (is(enaset, "ENAset")) {
    warning(paste0("Usage of ENAset objects will be deprecated ",
      "and potentially removed altogether in future versions."))

    enaset <- ena.set(enaset);
  }

  font.family = match.arg(font.family);

  plot = ENAplot$new(enaset,
                     title,
                     dimension.labels,
                     font.size,
                     font.color,
                     font.family,
                     scale.to = scale.to,
                     ...
                   );

  return(plot);
}

"""
ena_plot_group = """
##
#' @title Plot of ENA set groups
#'
#' @description Plot a point based on a summary statistic computed from a given method (typically, mean) for a set of points in a projected ENA space
#'
#' @details Plots a point based on a summary statistic for a group (typically, mean)
#'
#' @export
#'
#' @param enaplot \code{\link{ENAplot}} object to use for plotting
#' @param points A matrix or data.frame where columns contain coordinates of points in a projected ENA space
#' @param method A function for computing a summary statistic for each column of points
#' @param labels A character which will be the label for the group's point
#' @param colors A character, determines color of the group's point, default: enaplot$color
#' @param shape A character, determines shape of the group's point, choices:  square, triangle, diamond, circle, default: square
#' @param confidence.interval A character that determines how the confidence interval is displayed, choices: none, box, crosshair, default: none
#' @param outlier.interval A character that determines how outlier interval is displayed, choices: none, box, crosshair, default: none
#' @param label.offset character: top left (default), top center, top right, middle left, middle center, middle right, bottom left, bottom center, bottom right
#' @param label.font.size An integer which determines the font size for label, default: enaplot\$font.size
#' @param label.font.color A character which determines the color of label, default: enaplot\$font.color
#' @param label.font.family A character which determines font type, choices: Arial, Courier New, Times New Roman, default: enaplot\$font.family
#' @param show.legend Logical indicating whether to show the point labels in the in legend
#' @param legend.name Character indicating the name to show above the plot legend
#' @param ... Additional parameters
#'
#' @import magrittr
#'
#' @seealso \code{\link{ena.plot}}, \code{ena.plot.points}
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
#' unitNames = set$enadata$units
#'
#' ### Plot Condition 1 Group Mean
#' plot = ena.plot.group(plot, as.matrix(set$points$Condition$FirstGame), labels = "FirstGame",
#'     colors = "red", confidence.interval = "box")
#'
#' ### plot Condition 2 Group Mean
#' plot = ena.plot.group(plot, as.matrix(set$points$Condition$SecondGame), labels = "SecondGame",
#'     colors  = "blue", confidence.interval = "box")
#'
#' print(plot);
#'
#' @return The  \code{\link{ENAplot}} provided to the function, with its plot updated to include the new group point.
##
ena.plot.group <- function(
  enaplot,
  points = NULL,
  method = "mean",
  labels = NULL,
  colors = default.colors[1],
  shape = c("square", "triangle-up", "diamond", "circle"),
  confidence.interval = c("none", "crosshairs", "box"),
  outlier.interval = c("none", "crosshairs", "box"),
  label.offset = "bottom right",
  label.font.size = NULL,
  label.font.color = NULL,
  label.font.family = NULL,
  show.legend = T,
  legend.name = NULL,
  ...
) {
  shape = match.arg(shape);
  confidence.interval = match.arg(confidence.interval);
  outlier.interval = match.arg(outlier.interval);

  if(is.null(points)) {
    stop("Points must be provided.");
  } else if(is(points, "ena.points")) {
    points = remove_meta_data(points)
  }

  ### problem if outlier and confidence intervals selected for crosshair
  if(confidence.interval == "crosshairs" && outlier.interval == "crosshairs") {
    print("Confidence Interval and Outlier Interval cannot both be crosshair");
    print("Plotting Outlier Interval as box");
    outlier.interval = "box";
  }

  ### if group more than one row, combine to mean
  confidence.interval.values = NULL;
  outlier.interval.values = NULL;
  if(
    (is(points, "data.frame") || is(points, "matrix")) &&
    nrow(points) > 1
  ) {
    if(is.null(method) || method == "mean") {
      if(confidence.interval != "none") {
        confidence.interval.values = matrix(
          c(as.vector(t.test(points[,1], conf.level = 0.95)$conf.int), as.vector(t.test(points[,2], conf.level = 0.95)$conf.int)),
          ncol=2
        );
      }
      if(outlier.interval != "none") {
        outlier.interval.values = c(IQR(points[,1]), IQR(points[,2])) * 1.5;
      }

      if(length(unique(colors)) > 1) {
        points = t(sapply(unique(colors), function(color) colMeans(points[color == colors,]), simplify = T))
        colors = unique(colors)
        attr(enaplot, "means") <- length(attr(enaplot, "means")) + length(colors)
      } else {
        points = colMeans(points);
        attr(enaplot, "means") <- length(attr(enaplot, "means")) + 1
      }
    }
    else {
      if(confidence.interval != "none") warning("Confidence Intervals can only be used when method=`mean`")
      if(outlier.interval != "none") warning("Outlier Intervals can only be used when method=`mean`")

      points = apply(points, 2, function(x) do.call(method, list(x)) )
      attr(enaplot, "means") <- length(attr(enaplot, "means")) + 1
    }
  }

  enaplot %<>% ena.plot.points(
    points = points,
    labels = labels,
    colors = colors,
    shape = shape,
    confidence.interval = confidence.interval,
    confidence.interval.values = confidence.interval.values,
    outlier.interval = outlier.interval,
    outlier.interval.values = outlier.interval.values,
    label.offset = label.offset,
    label.font.size = label.font.size,
    label.font.color = label.font.color,
    label.font.family = label.font.family,
    show.legend = show.legend,
    legend.name = legend.name,
    ...
  )
  return(enaplot)

  #
  # group.layout = data.frame(dfDT.points);
  #
  # ### INTERVAL CALCULATIONS
  # error = NULL;
  # lines = list();
  #
  # if(confidence.interval == "crosshair") {
  #   ci.x = t.test(points.raw, conf.level = .95)$conf.int[1];
  #   ci.y = t.test(points.raw, conf.level = .95)$conf.int[2];
  #   error = list(
  #     x = list(type = "data", array = ci.x),
  #     y = list(type = "data", array = ci.y)
  #   )
  # } else if(outlier.interval == "crosshair") {
  #   oi.x = IQR(points.raw$V1) * 1.5;
  #   oi.y = IQR(points.raw$V2) * 1.5;
  #   error = list(
  #     x = list(type = "data", array = oi.x),
  #     y = list(type = "data", array = oi.y)
  #   )
  # }
  #
  # if(confidence.interval == "box") {
  #
  #   conf.ints = t.test(points.raw, conf.level = .95)$conf.int;
  #   dfDT.points[,c("ci.x", "ci.y") := .(conf.ints[1], conf.ints[2])]
  #
  #   #add cols for coordinates of CI lines
  #   dfDT.points[, c("ci.x1", "ci.x2", "ci.y1", "ci.y2") := .(V1 - ci.x, V1 + ci.x, V2 - ci.y, V2 + ci.y)]
  #
  #   lines.CI = apply(dfDT.points,1,function(x) {
  #     list(
  #       "type" = "square",
  #       "line" = list(
  #         width = 1,
  #         color = color,
  #         dash="dash"
  #       ),
  #       "xref" = "x",
  #       "yref" = "y",
  #       "x0" = x[['ci.x1']],
  #       "x1" = x[['ci.x2']],
  #       "y0" = x[['ci.y1']],
  #       "y1" = x[['ci.y2']]
  #     );
  #   });
  #   lines = lines.CI;
  # }
  # if(outlier.interval == "box") {
  #
  #   oi.x = IQR(points.raw$V1) * 1.5;
  #   oi.y = IQR(points.raw$V2) * 1.5;
  #
  #   dfDT.points[,c("oi.x", "oi.y") := .(oi.x, oi.y)]
  #
  #   #add cols for coordinates of CI lines
  #   dfDT.points[, c("oi.x1", "oi.x2", "oi.y1", "oi.y2") := .(V1 - oi.x, V1 + oi.x, V2 - oi.y, V2 + oi.y)]
  #
  #   lines.OI = apply(dfDT.points,1,function(x) {
  #     list(
  #       "type" = "square",
  #       "line" = list(
  #         width = 1,
  #         color = color,
  #         dash="dash"
  #       ),
  #       "xref" = "x",
  #       "yref" = "y",
  #       "x0" = x[['oi.x1']],
  #       "x1" = x[['oi.x2']],
  #       "y0" = x[['oi.y1']],
  #       "y1" = x[['oi.y2']]
  #     );
  #   });
  #
  #   lines = c(lines, lines.OI);
  # }
  #
  #
  # if(!is.null(error)) {
  #   #plot group w/ crosshair error bars
  #   enaplot$plot = plotly::add_trace(
  #     enaplot$plot,
  #     data = group.layout,
  #     type="scatter",
  #     x = ~V1, y = ~V2,
  #     mode="markers",
  #     marker = list(
  #       symbol =  shape,
  #       color = color,
  #       size = size
  #     ),
  #     error_x = error$x,
  #     error_y = error$y,
  #     showlegend = F,
  #     text = label,
  #     hoverinfo = "text+x+y"
  #   )
  # } else {
  #   #plot group w/o crosshair error bars
  #   enaplot$plot = plotly::add_trace(
  #     enaplot$plot,
  #     data = group.layout,
  #     type="scatter",
  #     x = ~V1, y = ~V2,
  #     mode="markers",
  #     marker = list(
  #       symbol =  shape,  #c(rep("circle",nrow(data)),rep("square", ifelse(!is.null(dfDT.groups), nrow(dfDT.groups), 0))),
  #       color = color,
  #       #size = c(rep(unit.size * unit.size.multiplier, nrow(data)), rep(group.size, ifelse(!is.null(dfDT.groups),nrow(dfDT.groups), 0)))
  #       size = size
  #     ),
  #     showlegend = F,
  #     text = label,
  #     hoverinfo = "text+x+y"
  #   )
  # }
  #
  # ##### WEIGHTING OFFSET
  # if(is.null(label.offset)) { label.offset = c(.05,.05) }
  # else label.offset = c(label.offset[1] * 0.1, label.offset[2] * 0.1)
  #
  # enaplot$plot = plotly::add_annotations(
  #   enaplot$plot,
  #   x = group.layout$V1[1] + label.offset[1],
  #   y = group.layout$V2[1] + label.offset[2],
  #   text = label,
  #   font = text.info,
  #   xref = "x",
  #   yref = "y",
  #   ax = label.offset[1],
  #   ay = label.offset[2],
  #   #xanchor = "left",
  #   showarrow = F
  # );
  #
  # enaplot$plot = plotly::layout(
  #   enaplot$plot,
  #   shapes = lines
  #   #annotations = label.info
  # )
  #
  # return(enaplot);
}
"""
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
  label.font.family = c("Arial", "Courier New"),
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
ena_rotate_by_mean = """
###
#' @title ENA Rotate by mean
#'
#' @description Computes a dimensional reduction from a matrix of points such that
#' the first dimension of the projected space passes through the means of two
#' groups in a the original space. Subsequent dimensions of the projected space
#' are computed using ena.svd
#'
#' @param enaset An \code{\link{ENAset}}
#' @param groups A list containing two logical vectors of length \code{nrow(ENA.set$ena.data$units)},
#' where each vector defines whether a unit is in one of the two groups whose means
#' are used to determine the dimensional reduction
#'
#' @export
#' @return \code{\link{ENARotationSet}}
###
ena.rotate.by.mean <- function(enaset, groups) {
  groups <- list(groups)
  groups <- groups[[1]]
  if (length(groups) < 1) {
    stop("Unable to rotate without 2 groups.")
  }

  if (!is(groups[[1]], "list")) {
    groups <- list(groups);
  }

  data <- as.matrix(enaset$line.weights)
  data <- scale(data, scale = F, center = T);

  col <- NULL
  vals <- NULL

  deflated.data <- data;
  i <- 1;
  weights <- matrix(0, nrow = ncol(deflated.data), ncol = length(groups))

  for (group in 1:length(groups)) {
    col <- group
    vals <- groups[[group]]

    col_one_vals <- deflated.data[vals[[1]], ]
    col_two_vals <- deflated.data[vals[[2]], ]
    col_one_means <- colMeans(as.matrix(col_one_vals))
    col_two_means <- colMeans(as.matrix(col_two_vals))
    col_mean_diff <- col_one_means - col_two_means

    col_mean_diff_sq <- col_mean_diff / sqrt(sum(col_mean_diff ^ 2))

    deflated.data <- deflated.data - (
                      deflated.data %*% col_mean_diff_sq
                    ) %*% t(col_mean_diff_sq)

    weights[, i] <- col_mean_diff_sq
    i <- i + 1;
  }
  defalted_data_svd <- orthogonal_svd(deflated.data, weights);

  colnames(defalted_data_svd) <- c(
    paste0("MR", as.character(1:length(groups))),
    paste0("SVD", as.character((length(groups) + 1):(ncol(defalted_data_svd))))
  )
  rownames(defalted_data_svd) <- colnames(as.matrix(enaset$line.weights))

  rotation_set <- ENARotationSet$new(
    node.positions = NULL,
    rotation = defalted_data_svd,
    codes = enaset$codes
  )
  return(rotation_set)
}

orthogonal_svd <- function(data, weights) {
  if (!is(data, "matrix")) {
    message("orthogonalSVD:  converting data to matrix")
    data <- as.matrix(data)
  }

  #Find the orthogonal transformation that includes W
  Q <- qr_ortho(weights)
  X.bar <- data %*% Q[, (ncol(weights) + 1):ncol(Q)]
  V <- prcomp(X.bar, scale. = F)$rotation
  if (is(V, "numeric")) {
    message("orthogonalSVD:  converting pca result to matrix")
    V <- matrix(V, nrow = length(V))
  }

  to_return <- (cbind(
    Q[, 1:ncol(weights)],
    Q[, (ncol(weights) + 1):ncol(Q)] %*% V
  ))

  return(to_return)
}

qr_ortho <- function(A) {
  return(qr.Q(qr(A), complete = T))
}

"""
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
ena_set_creator = """
#####
#'
#' @title Wrapper to generate an ENA model
#'
#' @description Generates an ENA model by constructing a dimensional reduction
#' of adjacency (co-occurrence) vectors as defined by the supplied
#' conversations, units, and codes.
#'
#' @details This function generates an ena.set object given a data.frame, units,
#' conversations, and codes. After accumulating the adjacency (co-occurrence)
#' vectors, computes a dimensional reduction (projection), and calculates node
#' positions in the projected ENA space. Returns location of the units in the
#' projected space, as well as locations for node positions, and normalized
#' adjacency (co-occurrence) vectors to construct network graphs. Includes options
#' for returning statistical tests between groups of units.
#'
#' @param data data.frame with containing metadata and coded columns
#' @param codes vector, numeric or character, of columns with codes
#' @param units vector, numeric or character, of columns representing units
#' @param conversation  vector, numeric or character, of columns to segment conversations by
#' @param metadata  vector, numeric or character, of columns with additional meta information for units
#' @param model character: EndPoint (default), AccumulatedTrajectory, SeparateTrajectory
#' @param weight.by "binary" is default, can supply a function to call (e.g. sum)
#' @param window MovingStanzaWindow (default) or Conversation
#' @param window.size.back Number of lines in the stanza window (default: 1)
#' @param include.meta [TBD]
#' @param groupVar vector, character, of column name containing group identifiers.
#' If column contains at least two unique values, will generate model using a means rotation (a dimensional reduction maximizing the variance between the means of the two groups)
#' @param groups vector, character, of values of groupVar column used for means rotation or statistical tests
#' @param runTest logical, TRUE will run a Student's t-Test and a Wilcoxon test for groups defined by the groups argument
#' @param ... Additional parameters passed to model generation
#'
#'
#' @return ena.set object
#####
ena.set.creator = function(
  data,
  codes,
  units,
  conversation,
  metadata = NULL,
  model = c("EndPoint", "AccumulatedTrajectory", "SeparateTrajectory"),
  weight.by = "binary",
  window = c("MovingStanzaWindow", "Conversation"),
  window.size.back = 1,
  # window.size.forward = 0,
  include.meta = TRUE,
  groupVar = NULL,
  groups = NULL,
  runTest = FALSE,
  # testType = c("nonparametric","parametric"),
  ...
) {
  data <- data.table::data.table(data)

  model = match.arg(model)
  window = match.arg(window)
  # testType = match.arg(testType)
  accum = ena.accumulate.data(
    units = data[, ..units, drop = FALSE],
    conversation = data[, ..conversation, drop = FALSE],
    metadata = data[, ..metadata, drop = FALSE],
    codes = data[, ..codes, drop = FALSE],
    window = window,
    window.size.back = window.size.back,
    # window.size.forward = window.size.forward,
    weight.by = weight.by,
    model = model,
    # mask = mask,
    include.meta = include.meta,
    ...
  );

  accum$model$raw.input <- as.data.table(data);
  accum$model$raw.input$ENA_UNIT <- merge_columns_c(accum$model$raw.input, units)
  group1 = NULL
  group2 = NULL
  group1.rows = NULL
  group2.rows = NULL

  set_params = list(...)
  set_params$enadata = accum

  ### make set if no group column is specified
  if(is.null(groupVar)) {
    if(runTest == TRUE) {
      warning("Group variable and groups not specified. Unable to run test")
    }
  }

  ### make set if group column is specified, but groups are not
  else if(is.null(groups) == TRUE) {
    unique.groups = unique(data[,groupVar])

    if(length(unique.groups) == 1) {
      warning("Group variable only contains one unique value. ENAset has been created without means rotation")

      if(runTest == TRUE) {
        warning("Multiple groups not specified. Unable to run test ")
      }
    }

    else{
      group1 = unique.groups[1]
      group2 = unique.groups[2]

      warning(paste0("No groups specified. Defaulting to means rotation using first two unique group values of group variable: ",group1," and ",group2))

      set_params$rotation.by = ena.rotate.by.mean
      set_params$rotation.params = list(accum$meta.data[[groupVar]] == group1, accum$meta.data[[groupVar]] == group2)

      if(runTest == TRUE) {
        warning(paste0("No groups specified. Running test on the first two unique group values of the group variable: ",group1," and ",group2))
      }
    }
  }
  else if(length(groups) == 1) {
    warning("Only one group value specified. ENAset has been created without means rotation")

    if(runTest == TRUE) {
      warning("Multiple groups not specified. Unable to run test")
    }
  }
  else {
    group1 = groups[1]
    group2 = groups[2]

    if(length(groups) > 2) {
      warning(paste0("Only two groups are allowed for means rotation. ENAset has been created using a means rotation on the first two groups given: ",group1," and ",group2))
    }

    groups.missing = groups[which(!groups %in% data[[groupVar]])]
    if(length(groups.missing) > 0) {
      stop(paste("Group column does not contain supplied group value(s): ", groups.missing))
    }

    if(runTest == TRUE) {
      if(length(groups) > 2) {
        warning(paste0("More than two groups specified. Running test on the first two groups: ",group1," and ",group2))
      }
    }
  }

  if(!any(is.null(c(group1, group2)))) {
    set_params$rotation.by = ena.rotate.by.mean
    set_params$rotation.params = list(accum$meta.data[[groupVar]] == group1, accum$meta.data[[groupVar]] == group2)

    group1.rows = accum$meta.data[[groupVar]] == group1
    group2.rows = accum$meta.data[[groupVar]] == group2
  }

  set = do.call(ena.make.set, set_params)

  if(
    runTest == TRUE &&
    !any(is.null(c(group1.rows, group2.rows)))
  ) {
    group1.dim1 = as.matrix(set$points)[group1.rows,1]
    group2.dim1 = as.matrix(set$points)[group2.rows,1]
    group1.dim2 = as.matrix(set$points)[group1.rows,2]
    group2.dim2 = as.matrix(set$points)[group2.rows,2]

    set$tests = list(
      wilcox.test = list(
        test.dim1 = wilcox.test(x = group1.dim1, y = group2.dim1),
        test.dim2 = wilcox.test(x = group1.dim2, y = group2.dim2)
      ),
      t.test = list(
        test.dim1 = t.test(x = group1.dim1, y = group2.dim1),
        test.dim2 = t.test(x = group1.dim2, y = group2.dim2)
      )
    )
  } else {
    set$tests = NULL
  }

  return(set)
}
"""
ena_svd = """
###
#' @title ENA SVD
#' @description ENA method computing a dimensional reduction of points in an ENA set using SVD
#' @param enaset An \code{\link{ENAset}}
#' @param ... Unused, necessary for ena.make.set
#' @export
###
ena.svd <- function(enaset, ...) {
  # to.norm = data.table::data.table(
  #   enaset$points.normed.centered,
  #   enaset$enadata$unit.names
  # )
  # to.norm = as.matrix(to.norm[,tail(.SD,n=1),.SDcols=colnames(to.norm)[which(colnames(to.norm) != "V2")],by=c("V2")][,2:ncol(to.norm)]);
  # pcaResults = pca_c(to.norm, dims = enaset$get("dimensions"));
  # pcaResults = pca_c(enaset$points.normed.centered, dims = enaset$get("dimensions"));

  # pts = enaset$model$points.for.projection[,!colnames(enaset$model$points.for.projection) %in% colnames(enaset$meta.data), with=F]
  pts = as.matrix(enaset$model$points.for.projection)
  pcaResults = prcomp(pts, retx=FALSE, scale=FALSE, center=FALSE, tol=0)

  ### used to be  enaset$data$centered$pca
  #enaset$rotation.set = pcaResults$pca;


  colnames(pcaResults$rotation) = c(
    paste('SVD',as.character(1:ncol(pcaResults$rotation)), sep='')
  );

  # rotationSet = ENARotationSet$new(rotation = pcaResults$pca, codes = enaset$codes, node.positions = NULL, eigenvalues = pcaResults$latent)
  rotationSet = ENARotationSet$new(
    rotation = pcaResults$rotation,
    codes = enaset$codes,
    node.positions = NULL,
    eigenvalues = pcaResults$sdev^2
  )
  return(rotationSet)
}

ena.svd.R6 <- function(enaset, ...) {
  pcaResults = prcomp(enaset$points.normed.centered, retx=FALSE,scale=FALSE,center=FALSE, tol=0)

  colnames(pcaResults$rotation) = c(
    paste('SVD',as.character(1:ncol(pcaResults$rotation)), sep='')
  );

  rotationSet = ENARotationSet$new(
    rotation = pcaResults$rotation, codes = enaset$codes,
    node.positions = NULL, eigenvalues = pcaResults$sdev^2
  )
  return(rotationSet)
}
"""
ena_unit_group = """
ena.unit.group = function(set,
  units, method = "mean", scale=T, name = "Group", scaleFactor = 1.0
  # ,keep.dimensions = 1:ncol(set$points.rotated)
) {
  runCIs <- function(pnts) {
    # pntRows = as.matrix(rep(T, nrow(set$points.rotated)))
    # if(length(units.by)>1) {
    #   pntRows = as.data.frame(set$enadata$units[[make.names(units.by[[1]])]]) == x;
    # }
    # pnts = as.matrix(set$points.rotated[pntRows,])
    # dim(pnts) = c(length(which(pntRows)),ncol(set$points.rotated))
    ci = matrix(NA, ncol=ncol(pnts),nrow=2)
    oi = rep(NA, ncol(pnts))
    if(nrow(pnts) > 1) {
      ci = t(apply(pnts,2,function(x) {
        tryCatch(t.test(x, conf.level = 0.95), error = function(e) list(conf.int = rep(x[1],2)))$conf.int
      })) * scaleFactor
      oi = apply(pnts, 2, function(x) { IQR(x) }) * 1.5 * scaleFactor
      # ci = t(matrix(c(
      #   tryCatch(t.test(pnts[, 1], conf.level = 0.95), error = function(e) list(conf.int = c(NA,NA)))$conf.int,
      #   tryCatch(t.test(pnts[, 2], conf.level = 0.95), error = function(e) list(conf.int = c(NA,NA)))$conf.int
      # ), nrow=2)) * scaleFactor;
      # oi = c(IQR(pnts[,1]), IQR(pnts[,2])) * 1.5 * scaleFactor;
    }
    list(
      ci = ci, #[keep.dimensions,],
      oi = oi #[keep.dimensions]
    )
  }

  runMean <- function(x) {
    group = list(
      "names" = as.vector(unique(x)),
      "points" = ena.group(set$points.rotated, x, method=mean),
      "line.weights" = ena.group(set$line.weights, x, method=method)
    )
    group$points = (as.vector(as.matrix(group$points)) * scaleFactor) #[keep.dimensions];

    colnames(group$points) <- NULL;
    if(method == "sum") {
      group$line.weights = group$line.weights * length(which(x == T));
    }
    # cis = runCIs(set$points.rotated[x,]);
    cis = runCIs(matrix(set$points.rotated[x,], ncol = ncol(set$points.rotated)));

    group$conf.ints = cis$ci;
    group$outlier.ints = cis$oi;
    group$line.weights = as.vector(as.matrix(group$line.weights))
    colnames(group$line.weights) = NULL;

    if(scale == T) {
      if(method == "sum") {
        group$line.weights = scales::rescale(group$line.weights, c(0.1,1));
      }
      group$edge.saturation = scales::rescale(group$line.weights, c(0.25,1));
      group$edge.opacity = scales::rescale(group$line.weights, c(0.3,1));
    }

    group
  }

  if(is.list(units)) {
    lapply(names(units), function(x) {
      unitsFnd = set$enadata$unit.names %in% unlist(units[x])
      grp = runMean(unitsFnd);
      grp$names = NULL;
      grp$name = unlist(as.character(x));
      grp$rle = list( lengths = length(which(unitsFnd == T)), values = x );
      grp;
    })
  } else {
    unitsFnd = set$enadata$unit.names %in% units
    grp = runMean(unitsFnd);
    grp$names = NULL;
    grp$name = name
    grp$rle = list( lengths = length(which(unitsFnd == T)), values = name );
    grp;
  }
}
"""
ena_unit_metadata = """
ena.unit.metadata = function(set) {

}
"""
ena_update_set = """
##
# @title Update ENA Set
#
# @description Update ENA set with a given ENA data object.
#
# @export
#
# @param ena.set \code{\link{ENAset}} from ena.make.set()
#
#
# @seealso \code{\link{ena.generate.set}}
#
# @examples
# \dontrun{
# #Given a ENA data object
# ena.update.set(\code{\link{ENAdata}})
# }
#
# @return \code{\link{ENAset}} class object
##
ena.update.set <- function(

) {

}
"""
ena_writeup = """
###
#' Calculate the correlations
#'
#' @description Calculate both Spearman and Pearson correlations for the
#' provided ENAset
#'
#' @param enaset ENAset to view methods of
#' @param tool c("rENA","webENA")
#' @param tool.version as.character(packageVersion(tool))
#' @param comparison character string representing the comparison used, c(NULL, "parametric", "non-parametric"). Default NULL
#' @param comparison.groups Groups that were used for the comparison
#' @param sig.dig Integer for the number of digits to round to
#' @param output_dir Where to save the output file
#' @param type c("file","stream") File will save to a file in output_dir, Stream returns the contents directly
#' @param theory Logical indicating whether to include theory in the writeup
#' @param methods Logical indicating whether to include methods in the writeup
#'
#' @export
#'
#' @return String representing the methods used to generate the model
ena.writeup <- function(
  enaset,
  tool = "rENA", tool.version = as.character(packageVersion(tool)),
  comparison = NULL, comparison.groups = NULL, sig.dig = 2,
  output_dir = getwd(), type = c("file","stream"), theory = T, methods = T
) {
  type = match.arg(type, choices = c("file","stream"), several.ok = FALSE)
  file = rmarkdown::render(system.file("rmd","methods.rmd", package="rENA"), output_dir = output_dir,
                    knit_root_dir = output_dir, intermediates_dir = output_dir, quiet = TRUE)
  if(type == "file") file
  else readChar(file, file.info(file)$size)
}

#' @title methods_report
#' @description Methods report for rmarkdwon
#' @param toc [TBD]
#' @param toc_depth [TBD]
#' @param fig_width [TBD]
#' @param fig_height [TBD]
#' @param keep_md [TBD]
#' @param md_extensions [TBD]
#' @param pandoc_args [TBD]
#'
#' @export
methods_report <- function(toc = FALSE,
                          toc_depth = 3,
                          fig_width = 5,
                          fig_height = 4,
                          keep_md = FALSE,
                          md_extensions = NULL,
                          pandoc_args = NULL) {

  # knitr options and hooks
  knitr <- rmarkdown::knitr_options(
    opts_chunk = list(dev = 'png',
                      dpi = 96,
                      fig.width = fig_width,
                      fig.height = fig_height)
  )

  # build pandoc args
  args <- c("--standalone")

  # table of contents
  args <- c(args, rmarkdown::pandoc_toc_args(toc, toc_depth))

  # pandoc args
  args <- c(args, pandoc_args)

  preserved_chunks <- character()

  # pre_processor <- function(metadata, input_file, runtime, knit_meta,
  #                           files_dir, output_dir) {
  #   preserved_chunks <<- extract_preserve_chunks(input_file, knitr::extract_raw_output)
  #   NULL
  # }

  # post_processor <- function(metadata, input_file, output_file, clean, verbose) {
  #   output_str <- readLines(output_file, encoding = 'UTF-8')
  #   output_res <- knitr::restore_raw_output(output_str, preserved_chunks)
  #   if (!identical(output_str, output_res))
  #     writeLines(enc2utf8(output_res), output_file, useBytes = TRUE)
  #   output_file
  # }

  # return output format
  rmarkdown::output_format(
    knitr = knitr,
    pandoc = rmarkdown::pandoc_options(to = "plain",
                            from = rmarkdown::from_rmarkdown(extensions = md_extensions),
                            args = args),
    keep_md = keep_md
    # ,pre_processor = pre_processor,
    # post_processor = post_processor
  )
}
"""
ENAdata = """####
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
      # BARRAcode{get()} - Return a read-only property
      # BARRApreformatted{  Example:
      #     get( x = 'file' )}
      # BARRApreformatted{  Parameters:
      #      x - Property to return. Defaults 'data', returns the original data
      ####
      get = function(x = "data") {
        return(private[[x]])
      },

      ####
      # BARRAcode{read()} - Return the accumulated data
      # BARRApreformatted{  Example:
      #     get( colnames = T, sep = " & " )}
      # BARRApreformatted{  Parameters:
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
        meta_cols_to_use <- meta_avail[apply(self$raw[, lapply(.SD, unique),
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
#' @field enaset - The BARRA code{BARRAlink{ENAset}} object from which the ENAplot was constructed
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
      # BARRAcode{get()} - Return a read-only property
      # BARRApreformatted{  Example:
      #     get( x = 'title' )}
      # BARRApreformatted{  Parameters:
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
ENArotation_set = """
#' ENARotationSet R6class
#
## Node positions
## Rotation matrix
#
#' @docType class
#' @importFrom R6 R6Class
#' @import data.table
#' @export
#
# @param windowSize Integer used to select the size of each stanza window within a conversation
# @param binary Logical, whether to convert code values to binary or allow for weigthed values
# @param unitsSelected deprecated
#
# @section Public ENARotationSet methods:
#######
ENARotationSet = R6::R6Class("ENARotationSet",
  public = list(

    #######
    ### Constructor - documented in main class declaration
    #######
    initialize = function(
      rotation,
      codes,
      node.positions,
      eigenvalues = NULL
    ) {
      self$node.positions = node.positions;
      self$rotation = rotation;
      self$codes = codes;
      if(!is.null(codes) && !is.null(self$node.positions)) {
       rownames(self$node.positions) = codes;
      }
      self$eigenvalues = eigenvalues;
    },

    ####
    ## Public Properties
    ####
      rotation = NULL,
      node.positions = NULL,
      codes = NULL,
      eigenvalues = NULL,
    ####
    ## END: Public Properties
    ####,

    ####
    ## Public Functions
    ####
      print = function(...) {
        args = list(...);
        fields = NULL;
        to.print = list();

        if (!is.null(args$fields)) {
          fields = args$fields
        } else if(!is.null(private$args$fields)) {
          fields = private$args$fields
        } else {
          fields = Filter(function(f) {
            cls = class(self[[f]]);
            !is(self[[f]], "function") && !is.null(self[[f]])
          }, names(get(class(self))$public_fields))
        }

        for(field in fields) {
          if(grepl("//", field)) {
            parts = Filter(function(f) { f!="" }, strsplit(field,"//")[[1]])
            to.print[[field]] = Reduce(function(o, i) { o[[i]] }, parts, self)
          } else {
            to.print[[field]] = self[[field]]
          }
        }
        return(to.print);
      }
    ####
    ## END: Public Functions
    ####
  ),
  private = list(
    #####
    ## Private Properties
    #####
      args = NULL
    #####
    ## END: Private Properties
    #####
  )
)
"""

ENAset = """####
#' ENAset R6class
#'
#' @docType class
#' @importFrom R6 R6Class
#' @import data.table
#' @export

#' @field enadata An BARRAcode{BARRAlink{ENAdata}} object originally used to create the set
#' @field points.raw A data frame containing accumulated adjacency (co-occurrence) vectors per unit
#' @field points.normed.centered A data frame of centered normed accumulated adjacency (co-occurrence) vectors for each unit
#' @field points.rotated A data frame of point positions for number of dimensions specified in ena.make.set (i.e., the centered, normed, and rotated data)
#' @field line.weights A data frame of connections strengths per unit (Data frame of normed accumu- lated adjacency (co-occurrence) vectors for each unit)
#' @field node.positions - A data frame of positions for each code
#' @field codes - A vector of code names
#' @field rotation.set - An BARRAcode{BARRAlink{ENARotationSet}} object
#' @field correlation - A data frame of spearman and pearson correlations for each dimension specified
#' @field variance - A vector of variance accounted for by each dimension specified
#' @field centroids - A matrix of the calculated centroid positions
#' @field function.call - The string representation of function called
#' @field function.params - A list of all parameters sent to function call
#'
####

ENAset = R6::R6Class("ENAset",
   public = list(
     #####
     ### Constructor - documented in main class declaration
     #####
     initialize = function(
       enadata,
       dimensions = 2,

       norm.by = fun_sphere_norm,

       rotation.by = ena.svd.R6,
       rotation.params = NULL,
       rotation.set = NULL,

       #center.data = center_data_c,    ### made local to run
       node.position.method = lws.positions.sq.R6,
       endpoints.only = T,
       ...
     ) {
       self$enadata <- enadata;

       private$dimensions <- dimensions;

       self$codes <- enadata$codes;

       self$function.call <- sys.call(-1);

       self$function.params$norm.by <- norm.by;    #was sphere_norm
       #self$function.params$center.data <- center.data;
       self$function.params$node.position.method <- node.position.method;    #was position.method
       self$function.params$rotation.by <- rotation.by;
       self$function.params$rotation.params <- rotation.params;
       self$function.params$rotation.set <- rotation.set;
       self$function.params$endpoints.only <- endpoints.only;

       private$args <- list(...);
     },
     #####
     ### END: Constructor
     #####

     #####
     ## Public Properties
     #####
     rotation_dists = NULL,  #leave for now - to be removed for a temp variable
     enadata = NULL,
     points.raw = NULL,    #was data$raw
     points.normed.centered = NULL,    #was data$centered$normed
     points.rotated = NULL,    #was data$centered$rotated
     points.rotated.scaled = NULL,
     points.rotated.non.zero = NULL,
     line.weights = NULL,   #was data$normed
     line.weights.non.zero = NULL,
     line.weights.unrotated = NULL,
     node.positions = NULL,  #was nodes$positions$scaled
     codes = NULL,
     rotation.set = NULL,   ## new - ENARotation object
     correlations = NULL,   #not formerly listed, comes from optimized node positions in egr.positions
     variance = NULL,     #was self$data$centered$latent
     centroids = NULL,
     function.call = NULL,     #new - string reping function call
     function.params = list(   #list containing parameters function was called with
       norm.by = NULL,
       node.position.method = NULL,
       rotation.by = NULL,
       rotation.params = NULL,
       endpoints.only = NULL
     ),
     #####
     ## END: Public Properties
     #####

     #####
     ## Public Functions
     #####

     ####
     # BARRAcode{process()} - Process the ENAset.
     # BARRApreformatted{}
     ####
     process = function() {
       return(private$run())
     },

     ####
     get.data = function(wh = c("normed","centered","rotated"), with.meta = T) {
       wh =  match.arg(wh);
       data = NULL;
       if( wh == "normed" ) {
         data = self$line.weights
       } else if ( wh == "centered" ) {
         data = self$points.normed.centered
       } else if ( wh == "rotated" ) {
         data = self$points.rotated
       }
       df.to.return = NULL;
       if(with.meta == T) {
         data.units = attr(data, opts$UNIT_NAMES);
         df.to.return = merge(
           data.table::data.table(
             data, data.units,
             ENA_UNIT=merge_columns_c(data.units, self$enadata$get("units.by")),
             TRAJ_UNIT=merge_columns_c(data.units, c(self$enadata$get("units.by"), self$enadata$get("trajectory.by")))
           ),
           self$enadata$add.metadata()
         )
       } else {
         df.to.return = data
       }
       df.to.return
     },

     ####
     # BARRAcode{get()} - Return a read-only property
     # BARRApreformatted{  Example:
     #     get( x = 'file' )}
     # BARRApreformatted{  Parameters:
     #      x - Property to return. Defaults to 'file', returning the original data}
     ####
     get = function(x = "enadata") {
       return(private[[x]])
     },
     print = function(...) {
       args = list(...);
       fields = NULL;
       to.print = list();

       if (!is.null(args$fields)) {
         fields = args$fields
       } else if(!is.null(private$args$fields)) {
         fields = private$args$fields
       } else {
         #fields = Filter(function(f) { (class(self[[f]]) != "function") }, names(get(class(self))$public_fields))
         fields = Filter(function(f) {
           cls = class(self[[f]]);
           !is(self[[f]], "function") && !is.null(self[[f]])
         }, names(get(class(self))$public_fields))
       }

       for(field in fields) {
         if(grepl("//", field)) {
           parts = Filter(function(f) { f!="" }, strsplit(field,"//")[[1]])
           to.print[[field]] = Reduce(function(o, i) { o[[i]] }, parts, self)
         } else {
           to.print[[field]] = self[[field]]
         }
       }
       return(to.print);
     }

     #####
     ## END: Public Functions
     #####
   ),

   private = list(
     #####
     ## Private Properties
     #####
     args = NULL,
     data.original = NULL,
     optim = NULL,

     #moved from public
     dimensions = 2,
     #####
     ## END: Private Properties
     #####

     #####
     ## Private Functions
     #####
     run = function() {       
       df = self$enadata$adjacency.vectors;

       ###
       # Backup of ENA data, this is not touched again.
       ###
       #private$data.original = df[,grep("adjacency.code", colnames(df)), with=F];
       private$data.original = df;
       ###
       # Copy of the original data, this is used for all
       # further operations. Unlike, `data.original`, this
       # is likely to be overwritten.
       ###
       self$points.raw = data.table::copy(private$data.original);

       ###
       # Normalize the raw data using self$function.params$norm.by,
       # which defaults to calling rENA::.sphere_norm
       ###
       self$line.weights = self$function.params$norm.by(self$points.raw);

       ###
       # Convert the string vector of code names to their corresponding
       # co-occurence names and set as colnames for the self$line.weights
       ##
       codeNames_tri = svector_to_ut(self$enadata$codes);

       colnames(self$line.weights) = codeNames_tri;
       # set the rownames to that of the original ENAdata file object
       rownames(self$line.weights) = rownames(df);

       attr(self$line.weights, opts$UNIT_NAMES) = attr(df, opts$UNIT_NAMES) #df[, .SD, with=T, .SDcols=self$enadata$get("unitsBy")];
       ###


       ###
       # Center the normed data
       # FIX - store as $data$centered
       ###
       #### ISSUE
       self$points.normed.centered = center_data_c(self$line.weights);

       colnames(self$points.normed.centered) = codeNames_tri;
       rownames(self$points.normed.centered) = rownames(df);
       attr(self$points.normed.centered, opts$UNIT_NAMES) = attr(self$enadata$adjacency.vectors.raw, opts$UNIT_NAMES)

       ###

       ###
       # Generate and Assign the rotation set
       ###
        if(!is.null(self$function.params$rotation.by) && is.null(self$function.params$rotation.set)) {
          self$rotation.set = do.call(self$function.params$rotation.by, list(self, self$function.params$rotation.params));
        } else if (!is.null(self$function.params$rotation.set)) {
          if(is(self$function.params$rotation.set, "ENARotationSet")) {
            print("Using custom rotation.set.")

            self$rotation.set = self$function.params$rotation.set;
          } else {
            stop("Supplied rotation.set is not an instance of ENARotationSet")
          }
        } else {
          stop("Unable to find or create a rotation set")
        }
       ###

       ###
       # Generated the rotated points
       ###
        self$points.rotated = self$points.normed.centered %*% self$rotation.set$rotation;
        private$dimensions = min(private$dimensions, ncol(self$points.rotated))
        attr(self$points.rotated, opts$UNIT_NAMES) = attr(self$points.normed.centered, opts$UNIT_NAMES);
       ###

       ###
       # Calculate node positions
       #  - The supplied methoed is responsible is expected to return a list
       #    with two keys, "node.positions" and "centroids"
       ###
        if(!is.null(self$rotation.set) && is.null(self$function.params$rotation.set)) {
          positions = self$function.params$node.position.method(self);
          if(all(names(positions) %in% c("node.positions","centroids"))) {
            self$node.positions = positions$node.positions
            self$centroids = positions$centroids

            self$rotation.set$node.positions = positions$node.positions
          } else {
            print("The node position method didn't return back the expected objects:")
            print("    Expected: c('node.positions','centroids')");
            print(paste("    Received: ",names(positions),sep=""));
          }
        } else if (!is.null(self$function.params$rotation.set)) {
          self$node.positions = self$function.params$rotation.set$node.positions
        } else {
          stop("Unable to determine the node positions either by calculating
                them using `node.position.method` or using a supplied
                `rotation.set`");
        }
       ###

       ###
       # Variance
       ###
       variance.of.rotated.data = var(self$points.rotated)
       diagonal.of.variance.of.rotated.data = as.vector(diag(variance.of.rotated.data))
       self$variance = diagonal.of.variance.of.rotated.data/sum(diagonal.of.variance.of.rotated.data)

       return(self);
     }
     #####
     ## END: Private Functions
     #####
   )
)
"""
namesToAdjacencyKey = """
##
#' @title Names to Adjacency Key
#'
#' @description Convert a vector of stringr, representing names of a square matrix, to an adjacency
#'
#' @details Returns a matrix of 2 rows by choose(length(vector), 2) columns
#'
#' @param vector Vector representing the names of a square matrix
#' @param upper_triangle Not Implemented
#'
#' @export
##
namesToAdjacencyKey <- function(vector, upper_triangle = TRUE) {
  upperTriIndices = triIndices(length(vector)) + 1;
  matrix(vector[upperTriIndices], nrow=2)
}
"""
ocpu_optimization = """
ocpu.optimization = function(
  set,
  maxit = 1000
) {

}
"""
RcppExports = """
# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' Merge data frame columns
#' @title Merge data frame columns
#' @description TBD
#' @param df Dataframe
#' @param cols Vector
#' @param sep Character seperator
#' @export
merge_columns_c <- function(df, cols, sep = ".") {
    .Call('_rENA_merge_columns_c', PACKAGE = 'rENA', df, cols, sep)
}

#' Upper Triangle from Vector
#'
#' @title vector to upper triangle
#' @description TBD
#' @param v [TBD]
#' @export
vector_to_ut <- function(v) {
    .Call('_rENA_vector_to_ut', PACKAGE = 'rENA', v)
}

svector_to_ut <- function(v) {
    .Call('_rENA_svector_to_ut', PACKAGE = 'rENA', v)
}

rows_to_co_occurrences <- function(df, binary = TRUE) {
    .Call('_rENA_rows_to_co_occurrences', PACKAGE = 'rENA', df, binary)
}

ref_window_df <- function(df, windowSize = 1, windowForward = 0, binary = TRUE, binaryStanzas = FALSE) {
    .Call('_rENA_ref_window_df', PACKAGE = 'rENA', df, windowSize, windowForward, binary, binaryStanzas)
}

ref_window_lag <- function(df, windowSize = 0L, binary = TRUE) {
    .Call('_rENA_ref_window_lag', PACKAGE = 'rENA', df, windowSize, binary)
}

#' Sphere norm
#' @title Sphere norm
#' @description TBD
#' @param dfM Dataframe
#' @export
fun_sphere_norm <- function(dfM) {
    .Call('_rENA_fun_sphere_norm', PACKAGE = 'rENA', dfM)
}

#' Non sphere norm
#'
#' @title Non sphere norm
#' @description TBD
#' @param dfM Dataframe
#' @export
fun_skip_sphere_norm <- function(dfM) {
    .Call('_rENA_fun_skip_sphere_norm', PACKAGE = 'rENA', dfM)
}

center_data_c <- function(values) {
    .Call('_rENA_center_data_c', PACKAGE = 'rENA', values)
}

triIndices <- function(len, row = -1L) {
    .Call('_rENA_triIndices', PACKAGE = 'rENA', len, row)
}

lws_lsq_positions <- function(adjMats, t, numDims) {
    .Call('_rENA_lws_lsq_positions', PACKAGE = 'rENA', adjMats, t, numDims)
}

# Register entry points for exported C++ functions
methods::setLoadAction(function(ns) {
    .Call('_rENA_RcppExport_registerCCallable', PACKAGE = 'rENA')
})

"""

rENA = """

#' @title rENA creates ENA sets
#' @description rENA is used to create and visualize network models of discourse and other phenomena from coded data using Epistemic Network Analysis (ENA). A more complete description of the methods will be provided with the next release. See also XXXXX
#' @name rENA
#' @importFrom Rcpp sourceCpp
#' @importFrom grDevices col2rgb
#' @importFrom grDevices hsv
#' @importFrom grDevices rgb2hsv
#' @importFrom methods is
#' @import stats
#' @import data.table
#' @import foreach
# @import plotly
#' @import utils
#' @import doParallel
#' @import parallel
# @import RcppRoll
# @import scales
# @import
# @import igraph
#' @useDynLib rENA
NULL

# @title Default rENA constants
# @description Default rENA constants
opts <- list (
  UNIT_NAMES = "ena.unit.names",
  TRAJ_TYPES = c("AccumulatedTrajectory", "SeparateTrajectory")
)

# @title Default colors used for plotting.
# @description Default colors for plotting
default.colors <- c(I("blue"), I("red"))


"""

RS_data = """

#' Coded Rescushell Chat Data
#'
#' A dataset containing sample chat data from the Rescushell Virtual
#' Internship
"RS.data"

"""

utils = """#' Find metadata columns
#'
#' @param x data.table (or frame) to search for columns of class ena.metadata
#'
#' @return logical vector
#' @export
find_meta_cols <- function(x) {
   sapply(x, is, class2 = "ena.metadata")
}

#' Find code columns
#'
#' @param x data.table (or frame) to search for columns of class ena.co.occurrence
#'
#' @return logical vector
#' @export
find_code_cols <- function(x) {
   grepl("adjacency.code", x = names(x)) | sapply(x, function(col) {
     is(col, class2 = "ena.co.occurrence")
   })
}

#' Find dimension columns
#'
#' @param x data.table (or frame) to search for columns of class ena.dimension
#'
#' @return logical vector
#' @export
find_dimension_cols <- function(x) {
   sapply(x, is, class2 = "ena.dimension")
}

#' Remove meta columns from data.table
#'
#' @param x [TBD]
#'
#' @return data.table withe columns of class ena.meta.data removed
#' @export
remove_meta_data <- function(x) {
   as.data.frame(x)[, !find_meta_cols(x), drop = F]
}

#' Extract metadata easily
#'
#' @param x [TBD]
#' @param i [TBD]
#'
#' @return [TBD]
#' @export
"$.ena.metadata" <- function(x, i) {
   #browser()
   parts <- unlist(strsplit(
               x = as.character(sys.call())[2], split = "//"
            ))[1:2]

   set <- get(parts[1], envir = parent.frame())
   m <- set[[parts[2]]][x == i, ]
   m
}

#' Extract line.weignts easily
#'
#' @param x [TBD]
#' @param i [TBD]
#'
#' @return [TBD]
#' @export
"$.line.weights" <- function (x, i) {
   vals <- x[[which(colnames(x) == i)]]

   vals
}

#' Extract points easily
#'
#' @param x [TBD]
#' @param i [TBD]
#'
#' @return [TBD]
#' @export
"$.ena.points" <- function (x, i) {
   vals <- x[[which(colnames(x) == i)]]

   vals
}
# "$.ena.plot" <- function(x, i) {
#  browser()
# }
# "[[.ena.plot" <- function(x, i) {
#  browser()
# }
#' @export
.DollarNames.ena.metadata <- function(x, pattern = "") {
   unique(x)
}

# "[.ena.matrix" = function(x, ...)
# {
#    browser()
#    original.class = class(x)[1]
#    class(x) = class(x)[-1]
#    x = x[...]
#
# #   y = as.data.frame(x)
# }

#' @export
summary.ena.set <- function(object, ...) {
   x <- object
   print_dims <- function(n = 2) {
      cat("     ", paste("Dimension", 1:n, collapse = "     "), " | ")
   }
   cat("Units: ", nrow(x$points), "          ")
   cat("Codes: ", length(x$rotation$codes), " | ")
    
   cat("Variance: | ")
   print_dims()
   cat("     ", paste(round(x$model$variance[1:2], 3), collapse = "     "), " |  |  ")

   cat("Eigenvalues: | ")
   print_dims()
   cat("      ", paste(round(
      x$rotation$eigenvalues[1:2], 3), collapse = "     "), "|  |   ")

   cat("Correlations: | ")
   cors <- ena.correlations(x)
   rownames(cors) <- paste("Dimension", 1:2)
   print(cors)
}
# as.data.frame.ena.connections <- function(x) {
#   class(x) = class(x)[-1]
#   y = as.data.frame(x)
#   y
# }
# format.co.occurrence = format.metadata = function(x, justify = "none") {
#   y = as.character(x)
#   format(y, justify = justify)
# }

#' Title
#'
#' @param x [TBD]
#' @param ... [TBD]
#' @param plot [TBD]
#' @param set [TBD]
#'
#' @return [TBD]
#' @export
print.ena.set <- function(x, ..., plot = FALSE, set = TRUE) {
   x.unclass <- unclass(x)

   if(
      !is.null(x.unclass$`_plot_op`) &&
      x.unclass$`_plot_op` == T
   ) {
      base::print(x.unclass$plots)
   }
   else {
      if(plot == FALSE) {
         x.unclass$plots <- NULL
      }
      base::print(x.unclass)
   }
}

#' Title
#'
#' @param x [TBD]
#' @param by [TBD]
#' @param model [TBD]
#' @param ... [TBD]
#'
#' @return [TBD]
#' @export
as_trajectory <- function(x,
   by = x$`_function.params`$conversation[1],
   model = c("AccumulatedTrajectory", "SeperateTrajectory"),
   ...
) {
   model = match.arg(model)
   orig_args = x$`_function.params`
   orig_args$model = model

   more_args <- list(...)
   for(arg in names(more_args)) {
      orig_args[[arg]] <- more_args[[arg]]
   }
   #c(mean, more.args[!names(more.args) %in% names(mean)])

   do.call(ena, orig_args)
}

#' Title
#'
#' @param x [TBD]
#' @param by [TBD]
#' @param ... [TBD]
#'
#' @return [TBD]
#' @export
project_in <- function(x, by = NULL, ...) {
   if(is.null(by)) {
      stop("A second parameter (ena.set or rotation.set) is required")
   }

   rotation.set <- NULL
   if(is(by, "ena.set")) {
      rotation.set <- by$rotation
   } else if(is(by, "ena.rotation.set")) {
      rotation.set <- by
   }

   if(!identical(x$rotation$adjacency.key, rotation.set$adjacency.key)) {
      stop("Rotation sets must have identical adjacency keys")
   }

   x$rotation.matrix <- rotation.set$rotation.matrix
   x$rotation$rotation.matrix <- rotation.set$rotation.matrix
   x$rotation$nodes <- rotation.set$nodes;
   x$rotation$eigenvalues <- rotation.set$eigenvalues

   points <- as.matrix(x$model$points.for.projection) %*% as.matrix(x$rotation.matrix)
   points.dt <- as.data.table(points)
   for (i in seq(ncol(points.dt))) {
    set(points.dt, j = i, value = as.ena.dimension(points.dt[[i]]))
   }
   if(grepl(x = x$model$model.type, pattern = "Trajectory")) {
    x$points <- cbind(x$trajectories, points.dt)
   } else {
    x$points <- cbind(x$meta.data, points.dt)
   }
   x$points <- as.ena.matrix(x$points, "ena.points")

   .return(x, invisible = T)
}

#' Title
#'
#' @param x [TBD]
#' @param on [TBD]
#'
#' @return [TBD]
#' @export
means_rotate <- function(x, on = NULL) {
   groupVar = NULL
   groups = NULL
   if(is.null(on)) {
      col_counts = as.numeric(x$model$raw.input[, lapply(.SD, function(s) {
                  length(unique(s))
               }),
               .SDcols = c(x$`_function.params`$units)
            ])
      groupVar = x$`_function.params`$units[order(col_counts) == 1]
      groups = levels(unique(x$model$raw.input[[groupVar]]))[1:2]
      # on_grps = list()
      # on_grps[[on]] = sapply(on_vals, function(v) {
      #    x$meta.data[[on]] == v
      # }, simplify = F)
   } else if(!is.null(names(on))) {
      groupVar = names(on)
      groups = on[[groupVar]]
   }

   if(is.null(groupVar) || is.null(groups)) {
      stop("Unable to determine groups for rotation.")
   }

   orig_args <- x$`_function.params`
   orig_args$groupVar = groupVar
   orig_args$groups = groups
   new_set <- do.call(ena, orig_args)
   new_set$plots <- x$plots
   invisible(new_set)
}

.return <- function(x, invisible = T, from_plot = F) {
   # browser()
   x$`_plot_op` = from_plot
# if() {
#       print(x$plots)
#    }

   if(invisible == T) {
      invisible(x)
   } else {
      return(x)
   }
}
"""

utils_classes = """

#' Re-class matrix as ena.matrix
#'
#' @param x data.frame, data.table, or matrix to extend
#' @param new.class Additional class to extend the matrix with, default: NULL
#'
#' @return Object of same st
#' @export
as.ena.matrix <- function(x, new.class = NULL) {
  class(x) = c(new.class, "ena.matrix", class(x))
  x
}

#' Re-class matrix as ena.metadata
#'
#' @param x data.frame, data.table, or matrix to extend
#'
#' @return Object of same st
#' @export
as.ena.metadata <- function(x) {
  if(is.factor(x)) {
    x = as.character(x)
  }
  class(x) = c("ena.metadata", "character") # This fails in the $.ena.metadata if is extending character, class(x))
  x
}
as.ena.code <- function(x) {
  if(is.factor(x)) {
    x = as.character(x)
  }
  class(x) = c("ena.code", class(x))
  x
}
as.ena.codes <- function(x) {
  if(is.factor(x)) {
    x = as.character(x)
  }
  class(x) = c("ena.codes", class(x))
  x
}
as.ena.co.occurrence <- function(x) {
  if(is.factor(x)) {
    x = as.character(x)
  }
  class(x) = c("ena.co.occurrence", class(x))
  x
}
as.ena.dimension <- function(x) {
  if(is.factor(x)) {
    x = as.character(x)
  }
  class(x) = c("ena.dimension", class(x))
  x
}
"""
utils_matrix = """
#' ENA line weights as matrix
#'
#' @param x ena.line.weights data.table to covert to matrix
#' @param ... additional arguments to be passed to or from methods
#' @param square [TBD]
#'
#' @return matrix
#' @export
as.matrix.ena.line.weights <- function(x, ..., square = FALSE) {
  args = list(...)

  # if(!is.null(args$square))
  #   square = args$square

  class(x) = class(x)[-1]
  rows = x[, !find_meta_cols(x), with = F]

  if(square) {
    upperTriSize = ncol(rows)
    number = ( (ceiling(sqrt(2*upperTriSize)) ^ 2) ) - (2*upperTriSize)
    codes = unique(unlist(sapply(colnames(rows), strsplit, split = " & ")))
    cm = sapply(seq(nrow(rows)), function(unit) {
      m = matrix(NA, number,  number, dimnames = list(codes, codes))
      m[upper.tri(m)] = as.numeric(rows[unit,])
      m
    }, simplify = F);
    return(cm)
  }
  else {
   as.matrix(remove_meta_data(rows), ...)
  }
}

#' ENA rotations as matrix
#'
#' @param x ena.rotation.matrix to conver to matrix
#' @param ... 	additional arguments to be passed to or from methods
#'
#' @return matrix
#' @export
as.matrix.ena.rotation.matrix <- function(x, ...) {
  class(x) = class(x)[-1]
  x = remove_meta_data(x)
  as.matrix(x, ...)
}

#' ENA points as matrix
#'
#' @param x ena.points to convert to a matrix
#' @param ... 	additional arguments to be passed to or from methods
#'
#' @return matrix
#' @export
as.matrix.ena.points <- function(x, ...) {
  class(x) = class(x)[-1]
  x = remove_meta_data(x)
  as.matrix(x, ...)
}

#' Matrix without metadata
#'
#' @param x Object to convert to  a matrix
#' @param ... 	additional arguments to be passed to or from methods
#'
#' @return matrix
#' @export
as.matrix.ena.matrix <- function(x, ...) {
  class(x) = class(x)[-1]
  x = remove_meta_data(x)
  as.matrix(x, ...)
}

#' ENA nodes as matrix
#'
#' @param x ena.nodes to convert to matrix
#' @param ... 	additional arguments to be passed to or from methods
#'
#' @return matrix
#' @export
as.matrix.ena.nodes <- function(x, ...) {
  class(x) = class(x)[-1]
  as.matrix(x[,-c("code")], ...)
}

#' ENA row connections as matrix
#'
#' @param x ena.row.connections to conver to a matrix
#' @param ... 	additional arguments to be passed to or from methods
#'
#' @return matrix
#' @export
as.matrix.row.connections <- function(x, ...) {
  class(x) = class(x)[-1]
  as.matrix(x[, sapply(x, is, class2="ena.co.occurrence"), with = F], ...)
}


#' ENA Connections as a matrix
#'
#' @param x ena.connections object
#' @param ... 	additional arguments to be passed to or from methods
# @param square Logical. If TRUE, each row is converted to a square matrix
# @param simplify Logical. If TRUE, returns back a single result as vector
# @param names Ignored
#'
#' @return If square is FALSE (default), a matrix with all metadata columns removed, otherwise a list with square matrices
#' @export
as.matrix.ena.connections <- function(x, ...) {
  class(x) = class(x)[-1]
  xx = remove_meta_data(x)
  rows = as.data.frame(x)[, !find_meta_cols(x), drop = F]

  args = list(...)
  if(is.null(args$square))
    square = F
  else
    square = args$square

  names = args$names

  if(is.null(args$simplify))
    simplify = ifelse(nrow(x) > 1, F, T)
  else
    simplify = args$simplify

  if(square) {
    upperTriSize = ncol(rows)
    number = ( (ceiling(sqrt(2*upperTriSize)) ^ 2) ) - (2*upperTriSize)
    codes = unique(unlist(sapply(colnames(rows), strsplit, split = " & ")))
    cm = sapply(seq(nrow(rows)), function(unit) {
      m = matrix(NA, number, number, dimnames = list(codes, codes))
      m[upper.tri(m)] = as.numeric(rows[unit,])
      m
    }, simplify = F)

    if(simplify) {
      cm = cm[[1]]
    } else {
      names(cm) = names
    }
  } else {
    cm = as.matrix(rows)
    rownames(cm) = names
  }

  cm
}
"""
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
zzz = """
.onLoad <- function(libname, pkgname) {
  globalVariables(c(".","ENA_ROW_IDX","ENA_UNIT","V1","V2","V3","ci.x","ci.y","e","handle","name","unit.groups","V","graph_from_data_frame","%>%","%<>%","X1","X2","dfDT.points","points.raw","lines","KEYCOL","ENA_CONV","..groupCol","..units","..metadata", "..codes", "..conversation"))
#   op <- options()
#   op.rENA <- list(
#     UNIT_NAMES = "ena.unit.names",
#     TRAJ_TYPES = c("accumulated","non-accumulated")
#   );
#
#   toset <- !("rENA" %in% names(op))
#   print(paste("ToSet:", toset));
#
#   if(toset) {
#     options(rENA = op.rENA)
#   }
#
#   invisible()
}
"""


string1 = lws_positions_sq + cohens_d + connection_matrix + ena + ena_accumulate_data + ena_accumulate_data_file + accumulate_data
string2 = ena_conversations + ena_correlations + ena_generate + ena_get_stats+ena_group+ena_make_set+ena_optimize_set+ena_plot+ena_plot_group+ena_plot_network
string3 = ena_plot_points+ena_plot_subtraction+ena_plot_trajectory+ena_plotter+ena_rotate_by_mean+ena_set+ena_set_creator
string4 = ena_svd+ena_unit_group+ena_unit_metadata+ena_update_set + ena_writeup
string5 = namesToAdjacencyKey+ocpu_optimization+RcppExports + RS_data + utils_classes + utils_matrix + utils_plot + zzz +rENA

string7 = ENArotation_set

string8 = ENAplot+ENAdata+utils+ENAset

stringfinal = string1+string2+string3+string4+string5+string8

stringr_c = STAP(stringfinal, "stringr_c")

print("Keys")

print(stringr_c._rpy2r.keys())



#print(type(df))

#stringr_c.ena_accumulate_data_file(RS_data, units_by = 1, conversations_by = 2, codes = "1")

#print(df)
#stringr_c.ena_make_set()

#stringr_c.ena_plot(stringr_c.ena_set)
