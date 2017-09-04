#' Function that allows to generate the standard form of the used code into the chosen model
#'
#' @importFrom R6 R6Class
#' @param  code the computational code (function of X and theta)
#' @param  X the matrix of forced variables
#' @param Yexp the vector of the experiments
#' @param model string of the model chosen ("model1","model2")
#' @return A function f(theta) is return
#' @examples
#' ### For the first model
#' X <- cbind(runif(3),runif(3))
#' code <- function(X,theta)
#' {
#'   return(X[,1]+theta*X[,2])
#' }
#' Yexp <- runif(3)
#' foo <- model(code,X,Yexp,"model1")
#' foo$fun(3,1)
#'
#' ### For the second model
#' X <- cbind(runif(5),runif(5),runif(5))
#' code <- function(X,theta)
#' {
#'   return(X[,1]+theta[,1]*X[,2]+theta[,2]^2*X[,3])
#' }
#' Yexp <- runif(5)
#' foo <- model(code,X,Yexp,"model2",opt.emul=list(p=2,n.emul=100,PCA=FALSE,binf=c(0,0),bsup=c(1,1)))
#' foo$fun(c(3,3),1)
#'
#' # with the PCA
#' foo <- model(code,X,Yexp,"model2",opt.emul=list(p=2,n.emul=100,PCA=TRUE,binf=c(0,0),bsup=c(1,1)))
#' foo$fun(c(3,3),1)
#' @export
model <- function(code,X,Yexp,model,opt.emul=list(p=1,n.emul=100,PCA=TRUE,binf=0,bsup=1))
{
  obj <- model.class$new(code,X,Yexp,model)
  switch(model,
         model1={
           obj = model1.class$new(code,X,Yexp,model)
           return(obj)
         },
         model2={
           obj = model2.class$new(code,X,Yexp,model,opt.emul)
           return(obj)
         }
  )
}




unscale.vector <- function(x,binf,bsup){
  # Function which unscale a vector between two bounds
  ###################### Inputs #######################
  # x the vector to unscale
  # binf the lower bound
  # bsup the upper bound
  ###################### Ouputs #######################
  # y the unscaled vector
  n = length(x)
  if (length(binf)!=length(bsup))
  {
    print("Please enter bounds of the same size")
  }
  if (length(binf)==1)
  {
    y=rep(1,n)
    for (i in 1:n){
      y[i]<- bsup*x[i]+(1-x[i])*binf
    }
  }else
  {
    y=rep(1,n)
    for (i in 1:n){
      y[i]<- bsup[i]*x[i]+(1-x[i])*binf[i]
    }
  }
  return(y)
}

unscale.matrix.diag <- function(M,binf,bsup){
  #### Funcion which unscale only the diagonal component of a matrix
  ################### Inputs #######################
  # M the matrix
  # binf the lower bounds
  # bsup the upper bounds
  ################### Outputs ######################
  # The normalized diagonal
  n <- dim(M)[1]
  T <- rep(1,n)
  for (i in 1:n){
    T[i] <- M[i,i]
  }
  T <- unscale.vector(T,binf,bsup)
  for (i in 1:n){
    M[i,i] <- T[i]
  }
  return(M)
}


unscale <- function(M,binf,bsup,diag=FALSE,sym=FALSE){
  ### Function which unscale un matrix or a vector
  #################### Inputs ####################
  # M the matrix or the vector
  # binf the lower bound
  # bsup the upper bound
  # option diag= False if we want to unscale the whole matrix
  # option sym=False if we do not have a symetric matrix
  if (diag==FALSE){
    if(sym==FALSE){
      if (is.matrix(M)==FALSE){
        if (length(binf)==length(M) & length(bsup)==length(M)){
          temp <- unscale.vector(M,binf,bsup)
        } else {
          temp <- unscale.vector(M,binf,bsup)
        }
      } else{
        N <- nrow(M)
        n <- ncol(M)
        if (n==1){
          temp <- unscale.vector(M,binf,bsup)
        } else {
          if (n!=length(binf) & N!=length(bsup)){
            print('The bound s sizes are not compatible. Please enter a upper bound and a lower bound of the same size than the matrix.')
          } else {
            temp <- matrix(nrow = N,ncol = n)
            for (i in 1:N){
              temp[i,] <- unscale.vector(M[i,],binf,bsup)
            }
          }
        }
      }
    } else {
      # Pour les matrice symétriques en cours...
    }
  } else {
    temp<-unscale.matrix.diag(M,binf,bsup)
  }
  return(temp)
}

