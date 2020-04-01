# -*- coding: UTF-8 -*-
import rpy2.robjects as robjects

from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import SignatureTranslatedAnonymousPackage as STAP
from rpy2.robjects import IntVector, Formula

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
      grep("V\\d+", colnames(dfDT.co.occurrences))
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
      grep("V\\d+", colnames(dfDT.summed.units))
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
ena_accumulate_data = """
##
#' @title Accumulate data from a data frame into a set of adjacency (co-occurrence) vectors
#'
#' @description This function initializes an ENAdata object, processing conversations from coded data to generate adjacency (co-occurrence) vectors
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
cohens_d = """
"""
cohens_d = """
"""
cohens_d = """
"""
cohens_d = """
"""
cohens_d = """
"""
cohens_d = """
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
default.colors <- c(I("blue"), I("red"))"

"""
utils ="""
#' Find metadata columns
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
               x = as.character(sys.call())[2], split = "\\$"
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
      cat("\t", paste("Dimension", 1:n, collapse = "\t"), "\n")
   }
   cat("Units: ", nrow(x$points), "\t\t")
   cat("Codes: ", length(x$rotation$codes), "\n")

   cat("Variance: \n")
   print_dims()
   cat("\t", paste(round(x$model$variance[1:2], 3), collapse = "\t\t"), "\n\n")

   cat("Eigenvalues: \n")
   print_dims()
   cat("\t", paste(round(
      x$rotation$eigenvalues[1:2], 3), collapse = "\t\t"), "\n\n")

   cat("Correlations: \n")
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






string = utils_plot + ena + zzz + ena_generate

stringr_c = STAP(string, "stringr_c")
print("Keys")
print(stringr_c._rpy2r.keys())

