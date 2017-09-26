#' A Reference Class to generate different calibration methods after generating a model from
#' model.class
#' @importFrom R6 R6Class
#'
#' @examples
#'
#' @export
estim.class <- R6Class(classname = "estim.class",
                 public = list(
                   LSE = function(theta,fun=fun)
                   {
                     Ytemp <- fun(theta[-length(theta)],theta[length(theta)])
                     return(sum((Ytemp-obj$Yexp)^2))
                   },
                   opt = function(fun)
                   {
                     return(optim(c(0,0),self$LSE,fun=fun)$par)
                   }
                 ))







