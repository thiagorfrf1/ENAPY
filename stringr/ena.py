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