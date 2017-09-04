#' A Reference Class to generate different calibration methods after generating a model from
#' model.class
#'
#' @include R6
#'
#' @examples
#' X <- cbind(runif(3),runif(3))
#'code <- function(X,theta)
#'{
#'  return(X[,1]+theta*X[,2])
#'}
#'Yexp <- runif(3)
#'test <- model(code,X,Yexp,"model1")
#'test2 <- estim.class$new()
#'test2$opt(test)
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


