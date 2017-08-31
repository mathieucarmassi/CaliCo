#' A Reference Class to generates differents model objects
#'
#' @description See the function blabla which produces an instance of this class
#'
#' This class comes with a set of methods, some of them being useful for the user:
#' See the documentation for blabla... Other methods
#'  should not be called as they are designed to be used during the optimization process.
#'
#' Fields should not be changed or manipulated by the user as they are updated internally
#' during the estimation process.
#' @field code a function which takes in entry X and theta
#' @field X the matrix of the forced variables
#' @field theta the vector of parameters to estimate
#' @field Yexp the experimental output
#' @field n the number of experiments
#' @field model the model choice see documentation
#' @seealso The function ....
model.class <- R6Class(classname = "gen",
                 public = list(
                   code     = NULL,
                   X        = NULL,
                   Yexp     = NULL,
                   n        = NULL,
                   model    = NULL,
                   initialize = function(code=NA,X=NA,Yexp=NA,model=NA)
                   {
                     self$code  <- code
                     self$X     <- X
                     self$Yexp  <- Yexp
                     self$n     <- length(Yexp)
                     self$model <- model
                     private$checkModels()
                   }
                 ))

model.class$set("private","checkModels",
        function()
        {
          if (self$model != "model1" & self$model != "model2")
          {
            stop('Please elect a correct model')
          }
        })

model.class$set("public","model1",
          function(theta,sig2)
            {
            return(self$code(self$X,theta)+rnorm(self$n,0,sqrt(sig2)))
          })

model.class$set("public","model2",
        function(theta,sig2,rho)
        {
          return(rho*self$code(self$X,theta)+rnorm(self$n,0,sqrt(sig2)))
        })

