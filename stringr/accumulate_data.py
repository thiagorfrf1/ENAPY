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