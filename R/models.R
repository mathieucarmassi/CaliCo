model <- R6Class(classname = "model",
                 public = list(
                   code  = NULL,
                   X     = NULL,
                   Yexp  = NULL,
                   n     = NULL,
                   initialize = function(code=NA,X=NA,Yexp=NA)
                   {
                     self$code  <- code
                     self$X     <- X
                     self$Yexp  <- Yexp
                     self$n     <- length(Yexp)
                   }
                 ))

model$set("public","model1",
          function(theta,sig2)
            {
            return(self$code(self$X,theta)+rnorm(self$n,0,sqrt(sig2)))
          })

# To do in estim (rename to calib). Make initialize function of the elected model and the data to infer the estimation.
# The user does not have to be to call the class model!!!

estim <- R6Class(classname = "estim",
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


X <- cbind(runif(3),runif(3))
code <- function(X,theta)
{
  return(X[,1]+theta*X[,2])
}
Yexp <- runif(3)

obj <- model$new(code,X,Yexp)
test <- estim$new()
test$LSE(c(1,1),obj=obj)
test$opt(obj)



