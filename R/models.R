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
gen <- R6Class(classname = "gen",
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
                     private$selection()
                   }
                 ))

gen$set("private","checkModels",
        function()
        {
          if (self$model != "model1" & self$model != "model2")
          {
            stop('Please elect a correct model')
          }
        })

gen$set("public","model1",
          function(theta,sig2)
            {
            return(self$code(self$X,theta)+rnorm(self$n,0,sqrt(sig2)))
          })

gen$set("public","model2",
        function(theta,sig2,rho)
        {
          return(rho*self$code(self$X,theta)+rnorm(self$n,0,sqrt(sig2)))
        })

gen$set("private","selection",
        function()
          {
          if (self$model== "model1")
          {
            funSelect <- self$model1
          } else
            if (self$model== "model2")
          {
            funSelect <- self$model2
          }
        })

# Example of generation a model
X <- cbind(runif(3),runif(3))
code <- function(X,theta)
{
  return(X[,1]+theta*X[,2])
}
Yexp <- runif(3)
test <- gen$new(code,X,Yexp,"model1")


# To do in estim (rename to calib). Make initialize function of the elected model and the data to infer the estimation.
# The user does not have to be to call the class model!!!

estim <- R6Class(classname = "estim",
                 inherit = gen,
                 public = list(
                   LSE = function(theta,obj=obj)
                   {
                     Ytemp <- obj$model1(theta[-length(theta)],theta[length(theta)])
                     return(sum((Ytemp-obj$Yexp)^2))
                   },
                   opt = function(obj)
                   {
                     return(optim(c(0,0),self$LSE,obj=obj)$par)
                   }
                 ))

# The class predict uses an object estim and new data to predict a new Y

calib <- function(code,X,Y)
{
  obj <- model$new(code,X,Y)
  est <- estim$new()
  return(est$opt(obj))
}

X <- cbind(runif(3),runif(3))
code <- function(X,theta)
{
  return(X[,1]+theta*X[,2])
}
Yexp <- runif(3)

calib(code,X,Yexp)

