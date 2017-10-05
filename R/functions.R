#' Generates \code{\link{model.class}} objects
#'
#' \code{model} is a function that allows us to generates a calibration model and its likelihood.
#'
#' There is four kind of models in calibration. They are properly defined in [1].
#'
#'
#' @param  code the computational code (function of X and theta)
#' @param  X the matrix of the forced variables
#' @param Yexp the vector of the experiments
#' @param model string of the model chosen ("model1","model2","model3","model4")
#' by default "model1" is choosen.
#' @return \code{model} returns a \code{model.class} object
#' @author M. Carmassi
#' @seealso \code{\link{model.class}},
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
#' foo$likelihood(3,1)
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
#' foo$likelihood(c(3,3),1)
#'
#' # with the PCA
#' foo <- model(code,X,Yexp,"model2",opt.emul=list(p=2,n.emul=100,PCA=TRUE,binf=c(0,0),bsup=c(1,1)))
#' foo$fun(c(3,3),1)
#'
#' ### For the third model
#' X <- cbind(runif(5),runif(5),runif(5))
#' code <- function(X,theta)
#' {
#'   return(X[,1]+theta[1]*X[,2]+theta[2]^2*X[,3])
#' }
#' Yexp <- runif(5)
#' foo <- model(code,X,Yexp,"model3")
#' foo$fun(c(3,3),c(0.2,0.2),1)
#' foo$likelihood(c(3,3),c(0.2,0.2),1)
#'
#'
#' ### For the fourth model
#' X <- cbind(runif(5),runif(5),runif(5))
#' code <- function(X,theta)
#' {
#'   return(X[,1]+theta[1]*X[,2]+theta[2]^2*X[,3])
#' }
#' Yexp <- runif(5)
#' foo <- model(code,X,Yexp,"model4",opt.emul=list(p=2,n.emul=100,PCA=FALSE,binf=c(0,0),bsup=c(1,1)))
#' foo$fun(c(3,3),c(0.2,0.2),1)
#' foo$likelihood(c(3,3),c(0.2,0.2),1)
#'
#' # with the PCA
#' foo <- model(code,X,Yexp,"model4",opt.emul=list(p=2,n.emul=100,PCA=TRUE,binf=c(0,0),bsup=c(1,1)))
#' foo$fun(c(3,3),c(0.2,0.2),1)
#' @export
model <- function(code,X,Yexp,model="model1",opt.emul=list(p=1,n.emul=100,PCA=TRUE,binf=0,bsup=1))
{
  library(R6)
  library(FactoMineR)
  library(DiceDesign)
  library(DiceKriging)
  switch(model,
         model1={
           obj = model1.class$new(code,X,Yexp,model)
           return(obj)
         },
         model2={
           obj = model2.class$new(code,X,Yexp,model,opt.emul)
           return(obj)
         },
         model3={
           obj = model3.class$new(code,X,Yexp,model)
           return(obj)
         },
         model4={
           obj = model4.class$new(code,X,Yexp,model,opt.emul)
           return(obj)
         }
  )
}


#' Generates \code{\link{prior.class}} objects
#'
#' \code{prior} is a function that allows us to generate one or a list of classes
#'  in which the charasteristics of each prior are defined
#'
#' The realized estimation is realized similarly as it is defined in [1]
#'
#'
#' @param  type.prior the vector of the prior types selected
#' @param opt.prior list of the hyperparameters relatives to the prior selected
#' @return \code{prior} returns a \code{prior.class} object
#' @author M. Carmassi
#' @seealso \code{\link{model.class}}, \code{\link{estim.class}}
#' @examples
#' ### Only one prior is wanted
#' ###### For a Gaussian Prior
#' foo <- prior(type.prior="gaussian",opt.prior=list(c(0.5,0.001)))
#' hist(foo$prior())
#'
#' ###### For a Gamma Prior
#' foo <- prior(type.prior="gamma",opt.prior=list(c(0.2,0.3)))
#' hist(foo$prior())
#'
#' ##### For an inverse-Gamma Prior
#' foo <- prior(type.prior="invGamma",opt.prior=list(c(0.2,0.3)))
#' hist(foo$prior())
#'
#' ### For several priors
#' foo <- prior(type.prior=c("gaussian","gamma"),opt.prior=list(c(0.5,0.001),c(0.2,0.3)))
#' hist(foo$Prior1$prior())
#'
#' @export
prior <- function(type.prior,opt.prior,log=FALSE)
{
  library(R6)
  n <- length(type.prior)
  if (n == 1)
  {
  switch(type.prior,
         gaussian = {
           obj = gaussian.class$new(type.prior,opt.prior,log)
           return(obj)
         },
         gamma={
           obj = gamma.class$new(type.prior,opt.prior,log)
           return(obj)
         },
         invGamma={
           obj = invGamma.class$new(type.prior,opt.prior,log)
           return(obj)
         },
         unif={
           obj = unif.class$new(type.prior,opt.prior,log)
           return(obj)
         }
  )
  } else
  {
    NAmes <- c("Prior1")
    res <- list()
    for (i in 1:n)
    {
      if (i>1){NAmes <- cbind(NAmes,paste("Prior",i,sep=""))}
      switch(type.prior[i],
             gaussian = {
               obj = gaussian.class$new(type.prior[i],opt.prior[[i]],log)
             },
             gamma={
               obj = gamma.class$new(type.prior[i],opt.prior[[i]],log)
             },
             invGamma={
               obj = invGamma.class$new(type.prior[i],opt.prior[[i]],log)
             },
             unif={
               obj = unif.class$new(type.prior[i],opt.prior[[i]],log)
             }
      )
      res[[i]] <- obj
      names(res) <- NAmes
    }
    return(res)
  }
}


#' Generates \code{\link{estim.class}} objects
#'
#' \code{estim} is a function that allows us to generate a class in which the estimation is
#' done
#'
#' The realized estimation is realized similarly as it is defined in [1]
#'
#' @useDynLib calibrationCode
#' @importFrom Rcpp evalCpp
#'
#' @param  md a \code{\link{model.class}} object
#' @param pr a \code{\link{prior.class}} object
#' @param x data for calibration
#' @return opt list of options for the inference
#' @author M. Carmassi
#' @seealso \code{\link{model.class}}, \code{\link{prior.class}}
#' @examples
#' ### For the first model
#' X <- cbind(runif(3),runif(3))
#' code <- function(X,theta)
#' {
#'   return(X[,1]+as.vector(theta)*X[,2])
#' }
#' Yexp <- runif(3)
#' test <- estim(code,X,Yexp,model="model1",type.prior=c("gaussian","gamma"),log=TRUE,opt.prior=list(c(3,1),c(0.2,0.3)),opt.estim=list(Ngibbs=3000,Nmh=10000,thetaInit=c(3,1),k=c(0.1,0.1),sig=diag(2)))
#'
#' @export
estim <- function(code,X,Yexp,model="model1",type.prior,log=FALSE
                  ,opt.emul=list(p=1,n.emul=100,PCA=TRUE,binf=0,bsup=1),opt.prior,opt.estim)
{
  md <<- model(code,X,Yexp,model,opt.emul)
  pr <- prior(type.prior,opt.prior,log=TRUE)
  binf <- pr[[1]]$binf
  bsup <- pr[[1]]$bsup
  for (i in 2:length(type.prior))
  {
    binf <- cbind(binf,pr[[i]]$binf)
    bsup <- cbind(bsup,pr[[i]]$bsup)
  }
  if (length(type.prior) == 1)
  {
    logTest <- function(theta,sig2){return(log(md$likelihood(theta,sig2))+pr$prior(theta))}
  } else
  {
    logTest <- function(theta,sig2)
    {
      s <- 0
      for (i in 1:(length(theta)))
      {
        s <- s + pr[[i]]$prior(theta[i])
      }
      s <- s + pr[[(length(theta)+1)]]$prior(sig2)
      return(log(md$likelihood(theta,sig2)) + s)
    }
  }
  res <- MetropolisHastingsCpp(md$fun,opt.estim$Ngibbs,opt.estim$Nmh,opt.estim$thetaInit,
                               opt.estim$k,opt.estim$sig,Yexp,binf,bsup,logTest)
  return(res)
}


#' Function which unscale a vector between two bounds
#'
#' @param  x the vector to unscale
#' @param  binf the lower bound
#' @param bsup the upper bound
#' @return y the vector unscaled
#' @examples
#' X <- runif(3)
#' Y <-unscale.vector(X,c(10,10,1O),c(15,15,15))
#' print(Y)
#' @export
unscale.vector <- function(x,binf,bsup){
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


#' Funcion which unscale only the diagonal component of a matrix
#'
#' @param  M the matrix
#' @param  binf the lower bound
#' @param bsup the upper bound
#' @return the normalized diagonal
#' @examples
#' X <- ones(3)*runif(3)
#' Y <-unscale.matrix.diag(X,c(10,10,1O),c(15,15,15))
#' print(Y)
#' @export
unscale.matrix.diag <- function(M,binf,bsup){
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


#' Function which unscale un matrix or a vector
#'
#' @param M the matrix or the vector
#' @param  binf the lower bound
#' @param bsup the upper bound
#' @param diag default value False if we want to unscale the whole matrix
#' @param sym default value False if we do not have a symetric matrix
#' @return the unscaled vector or matrix
#' @examples
#' X <- ones(3)*runif(3)
#' Y <- unscale(X,c(10,10,1O),c(15,15,15))
#' print(Y)
#' @export
unscale <- function(M,binf,bsup,diag=FALSE,sym=FALSE){
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
            print('The bound s sizes are not compatible.
                  Please enter a upper bound and a lower bound of the same size than the matrix.')
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

