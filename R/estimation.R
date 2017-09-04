estim.class <- R6Class(classname = "estim.class",
                 inherit = model.class,
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


