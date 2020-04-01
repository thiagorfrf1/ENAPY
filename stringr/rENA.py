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