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
                             x.new     = NULL,
                             lowerPred = NULL,
                             upperPred = NULL,
                             meanPred  = NULL,
                             initialize = function(modelfit,x.new)
                             {
                               self$modelfit    <- modelfit
                               self$x.new       <- x.new
                               self$checkNewdata(x.new)
                               self$predictCal()
                             },
                             predictCal = function()
                             {
                               Dim <- length(self$modelfit$type.prior)
                               if (is.matrix(self$modelfit$X)==TRUE)
                                 {
                                  l <- nrow(self$x.new)
                                 } else
                                 {
                                  l <- length(self$x.new)
                                 }
                               MAP <- self$modelfit$MAP
                               COV <- cov(self$modelfit$out$THETA[-c(1:self$modelfit$burnIn),])
                               s   <- mvrnorm(100,MAP,COV)
                               res <- matrix(nr=l,nc=100)
                               for (i in 1:100)
                               {
                                 res[,i] <- self$modelfit$md$pred(s[i,1:(Dim-1)],s[i,Dim],x.new)$y
                               }
                               self$lowerPred <- apply(res,2,quantile,probs=0.05)
                               self$upperPred <- apply(res,2,quantile,probs=0.95)
                               self$meanPred  <- apply(res,2,quantile,probs=0.5)
                               self$print()
                             },
                             print = function()
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
                             plot = function(rdata=NULL)
                             {
                               ggdata <- data.frame(lower=self$lowerPred,upper=self$upperPred,x=x.new,y=self$meanPred)
                               p <- ggplot(ggdata)+geom_line(aes(y=y,x=x))+
                                 geom_ribbon(aes(ymin=lower, ymax=upper, x=x), alpha = 0.3)
                               return(p)
                              },
                             checkNewdata = function(x.new)
                             {
                               if (is.matrix(x.new)==TRUE)
                               {
                                 if (ncol(x.new)!=ncol(self$modelfit$X))
                                 {
                                   stop("Enter the same dimension for X")
                                 }
                               } else
                               {
                                 if (length(x.new)!= length(self$modelfit$X))
                                 {
                                   stop("Enter the same dimension for X")
                                 }
                               }
                             }
                            ))



