#' A Reference Class to generate different calibration methods after generating a model from
#' model.class
#'
#' @examples
#'
#' @export
estim.class <- R6::R6Class(classname = "estim.class",
                 public = list(
                    logTest   = NULL,
                    opt.estim = NULL,
                    initialize = function(logTest=NA,opt.estim=NA)
                    {
                      self$logTest   <- logTest
                      self$opt.estim <- opt.optim

                    }
                   ))

estim.model1.class <- R6::R6Class(classname = "estim.model1.class",
                                  inherit = estim.class,
                                  public = list(

                                  ))





