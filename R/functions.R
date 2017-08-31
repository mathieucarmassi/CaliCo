#' Function that allows to generate the standard form of the used code into the chosen model
#'
#' @param  code the computational code (function of X and theta)
#' @param  X the matrix of forced variables
#' @param Yexp the vector of the experiments
#' @param model string of the model chosen ("model1","model2")
#' @return A function f(theta) is return
#' @examples
#' X <- cbind(runif(3),runif(3))
#'code <- function(X,theta)
#'{
#'  return(X[,1]+theta*X[,2])
#'}
#'Yexp <- runif(3)
#'test <- model(code,X,Yexp,"model1")
#'test(2,3)
#' @export
model <- function(code,X,Yexp,model)
{
  obj <- gen$new(code,X,Yexp,model)
  if (model == "model1")
  {
    return(obj$model1)
  } else
    if (model == "model2")
    {
      return(obj$model2)
    }
}
