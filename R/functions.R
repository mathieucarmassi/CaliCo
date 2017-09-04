#' Function that allows to generate the standard form of the used code into the chosen model
#'
#' @param  code the computational code (function of X and theta)
#' @param  X the matrix of forced variables
#' @param Yexp the vector of the experiments
#' @param model string of the model chosen ("model1","model2")
#' @return A function f(theta) is return
#' @examples
#' X <- cbind(runif(3),runif(3))
#' code <- function(X,theta)
#' {
#'   return(X[,1]+theta*X[,2])
#' }
#' Yexp <- runif(3)
#' foo <- model(code,X,Yexp,"model1")
#' foo$fun(3,1)
#' @export
model <- function(code,X,Yexp,model,n.emul=100)
{
  obj <- model.class$new(code,X,Yexp,model)
  switch(model,
         model1={
           obj = model1.class$new(code,X,Yexp,model)
           return(obj)
         },
         model2={
           obj = model2.class$new(code,X,Yexp,model,n.emul)
           return(obj)
         }
  )
}

