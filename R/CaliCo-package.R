#' Bayesian calibration for computational codes
#'
#' The CaliCo package provides three categories of important functions:
#' \code{\link{model}}, \code{\link{prior}}, \code{\link{calibrate}} and \code{\link{prediction}}.
#'
#' @useDynLib CaliCo
#' @importFrom R6 R6Class
#' @importFrom stats rnorm
#' @import ggplot2 DiceKriging DiceDesign FactoMineR coda parallel testthat MASS
#'
#' @details
#' Package: CaliCo
#'
#' Type:    Package
#'
#' Version: 0.1.0
#'
#' Date:    2017-11-07
#'
#' License: GPL-2 | GPL-3
#'
#' @docType package
#' @author Mathieu Carmassi
#' @author Maintainer: \email{mathieu.carmassi@gmail.com}
#' @references Bachoc et al. (2014) <arXiv:1301.4114v2>
#' @references Bayarri et al. (2007 b) <doi:10.1198/004017007000000092>
#' @references Carmassi et al. (2018) <arXiv:1801.01810>
#' @references Cox et al. (2001) <doi:10.1016/S0167-9473(00)00057-8>
#' @references Hastings (1970) <doi:10.2307/2334940>
#' @references Higdon et al. (2004) <doi:10.1137/S1064827503426693>
#' @references Kennedy et al. (2001) <doi:10.1111/1467-9868.00294>
#' @references Kennedy et al. (2001b) <doi:10.1.1.28.2835>
#' @references Liu et al. (2009) <doi:10.1214/09-BA404>
#' @references Roustant et al. (2012) <hal:00495766v1>
#' @references Sacks et al. (1989) <doi:10.1214/ss/1177012413>
#' @references Santner et al. (2003, ISBN:978-1-4757-3799-8)
#' @examples
#' # Introduction to CaliCo
#' \dontrun{vignette("CaliCo-introduction")}
#' @name CaliCo
NULL
