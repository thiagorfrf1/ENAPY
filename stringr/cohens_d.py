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