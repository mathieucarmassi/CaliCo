#' Function that realize the prediction after calibration
#'
#'
#'
#'
#'
#' @export
prediction.class <- R6::R6Class(classname = "predict.class",
                           public = list(
                             modelfit = NULL,
                             x.new    = NULL,
                             initialize = function(modelfit,x.new)
                             {
                               self$modelfit <- modelfit
                               self$x.new    <- x.new
                               self$predictCal()
                             },
                             predictCal = function()
                             {
                               Dim <- length(self$modelfit$type.prior)
                               MAP <- self$modelfit$MAP
                               COV <- cov(self$modelfit$out$THETA[-c(1:self$modelfit$burnIn),])
                               s   <- mvrnorm(100,MAP,COV)
                               res <- matrix(nr=100,nc=1)
                               for (i in 1:100)
                               {
                                 res[i] <- self$modelfit$md$fun(MAP[i,1:(Dim-1)],MAP[i,Dim])$y
                               }
                               lowerPred <- apply(res,2,quantile,probs=0.05)
                               upperPred <- apply(res,2,quantile,probs=0.95)
                               meanPred  <- apply(res,2,quantile,probs=0.5)
                               self$print()
                             },
                             print = function()
                               function()
                               {
                                 cat("Call:\n")
                                 print(self$model)
                                 cat("\n")
                                 cat("With the function:\n")
                                 print(self$code)
                                 cat("\n")
                                 cat("Converged")
                                 cat("\n\n")
                                 cat("MAP estimator:\n")
                                 print(round(self$modelfit$MAP,5))
                                 cat("\n")
                                 cat("Mean:\n")
                                 print(round(self$modelfit$Mean,5))
                               },
                             plot = function()
                             {
                                return(0)
                             }
                            ))



