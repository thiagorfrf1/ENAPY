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