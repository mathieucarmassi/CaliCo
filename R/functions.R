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
#' X <- seq(0,1,length.out=100)
#' code <- function(X,theta)
#' {
#'   return((6*X-2)^2*sin(theta*X-4))
#' }
#' Yexp <- code(X,11)+rnorm(100,0,0.1)
#' # Generate the model
#' foo <- model(code,X,Yexp,"model1")
#' # Plot the results
#' foo$plot(11,0.1,X)
#' # Print the results for new data X
#' foo$fun(11,0.1)
#' # Get acces to the likelihood
#' foo$likelihood(11,0.1)
#'
#' ### For the second model
#' X <- seq(0,1,length.out=25)
#' code <- function(X,theta)
#' {
#'   return((6*X-2)^2*sin(theta*X-4))
#' }
#' Yexp <- code(X,11)+rnorm(25,0,0.1)
#' # Generate the model with setup for the Gaussian Process
#' foo <- model(code,X,Yexp,"model2",opt.emul=list(p=1,n.emul=50,PCA=FALSE),binf=8,bsup=14)
#' # Plot the model
#' foo$plot(11,0.1,X,points=TRUE)
#'
#' # with the PCA in stand by
#' # foo <- model(code,X,Yexp,"model2",opt.emul=list(p=2,n.emul=100,PCA=TRUE),binf=c(0,0),bsup=c(1,1))
#' # foo$fun(c(3,3),1)
#'
#' ### For the third model
#' X <- seq(0,1,length.out=100)
#' code <- function(X,theta)
#' {
#'   return((6*X-2)^2*sin(theta*X-4))
#' }
#'
#' Yexp <- code(X,11)+rnorm(100,0,0.1)
#' foo <- model(code,X,Yexp,"model3")
#' foo$plot(11,c(50,1),0.1,X)
#' foo$likelihood(11,c(50,1),0.1)
#'
#' ### For the fourth model
#' X <- seq(0,1,length.out=100)
#' code <- function(X,theta)
#' {
#'   return((6*X-2)^2*sin(theta*X-4))
#' }
#' Yexp <- code(X,11)+rnorm(100,0,0.1)
#' foo <- model(code,X,Yexp,"model4",opt.emul=list(p=1,n.emul=60,PCA=FALSE),binf=8,bsup=14)
#' foo$plot(11,c(50,1),0.1,X,points=FALSE)
#'
#' foo$likelihood(11,c(50,1),0.1)
#'
#' # with the PCA in stand by
#' # foo <- model(code,X,Yexp,"model4",opt.emul=list(p=2,n.emul=100,PCA=TRUE),binf=c(0,0),bsup=c(1,1))
#' # foo$fun(c(3,3),c(0.2,0.2),1)
#' @export
model <- function(code,X,Yexp,model="model1",opt.emul=list(p=1,n.emul=100,PCA=TRUE),binf=0,bsup=1)
{
  library(R6)
  library(FactoMineR)
  library(DiceDesign)
  library(DiceKriging)
  switch(model,
         model1={
           obj = model1.class$new(code,X,Yexp,model,binf,bsup)
           return(obj)
         },
         model2={
           obj = model2.class$new(code,X,Yexp,model,opt.emul,binf,bsup)
           return(obj)
         },
         model3={
           obj = model3.class$new(code,X,Yexp,model,binf,bsup)
           return(obj)
         },
         model4={
           obj = model4.class$new(code,X,Yexp,model,opt.emul,binf,bsup)
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
#' @import ggplot2
#' @import gridExtra
#' @param  type.prior the vector of the prior types selected
#' @param opt.prior list of the hyperparameters relatives to the prior selected
#' @return \code{prior} returns a \code{prior.class} object
#' @author M. Carmassi
#' @seealso \code{\link{model.class}}, \code{\link{estim.class}}
#' @examples
#' ### Only one prior is wanted
#' ###### For a Gaussian Prior
#' foo <- prior(type.prior="gaussian",opt.prior=list(c(0.5,0.001)))
#' foo$plot()
#'
#' ###### For a Uniform Prior
#' foo <- prior(type.prior="unif",opt.prior=list(c(0,1)))
#' foo$plot()
#'
#' ###### For a Gamma Prior
#' foo <- prior(type.prior="gamma",opt.prior=list(c(5,1)))
#' foo$plot()
#'
#' ##### For an inverse-Gamma Prior
#' foo <- prior(type.prior="invGamma",opt.prior=list(c(0.2,0.3)))
#' foo$plot()
#'
#' ### For several priors
#' foo <- prior(type.prior=c("gaussian","gamma"),opt.prior=list(c(0.5,0.001),c(5,1)))
#' grid.arrange(foo$Prior1$plot(),foo$Prior2$plot(),nrow=2)
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
#' ####### For the first model
#' ### The data set
#' X <- seq(0,1,length.out=100)
#' ### The code to calibrate
#' code <- function(X,theta)
#' {
#'   return((6*X-2)^2*sin(theta*X-4))
#' }
#' ### Simulated data
#' Yr   <- code(X,11)
#' Yexp <- Yr+rnorm(100,0,0.1)
#' ### Definition of the nature of the priors
#' type.prior=c("gaussian","gamma")
#' ### Definition of the prior hyperparameters
#' opt.prior=list(c(11,3),c(2,0.1))
#' ### Definition of the estimation option
#' opt.estim=list(Ngibbs=400,Nmh=1000,thetaInit=c(11,0.1),k=c(5e-4,5e-4),sig=diag(2))
#' ### Definition of the emulation options (for Model2 and Model4 exclusively)
#' opt.emul=list(p=1,n.emul=150,PCA=FALSE)
#'
#' test <- estim(code,X,Yr,Yexp,model="model1",type.prior,log=TRUE,opt.prior,opt.estim)
#' test2 <- estim(code,X,Yr,Yexp,model="model2",type.prior,log=TRUE,opt.prior,opt.estim,opt.emul)
#' test3 <- estim(code,X,Yr,Yexp,model="model3",type.prior,log=TRUE,opt.prior,opt.estim)
#' test4 <- estim(code,X,Yr,Yexp,model="model4",type.prior,log=TRUE,opt.prior,opt.estim,opt.emul)
#' test$plot()
#' test2$plot()
#' test3$plot()
#' test4$plot()
#'
#' # With leave one out cross validation
#' opt.valid =list(n.CV=10,k=NULL)
#' test <- estim(code,X,Yr,Yexp,model="model1",type.prior,log=TRUE,opt.prior,opt.estim,type.valid="loo",opt.valid=opt.valid)

#'
#'
#' #### with two parameters
#' X <- seq(0,1,length.out=100)
#' code <- function(X,theta)
#' {
#'   return((theta[2]*X-2)^2*sin(theta[1]*X-4))
#' }
#' ### Simulated data
#' Yr   <- code(X,c(11,6))
#' Yexp <- Yr+rnorm(100,0,0.1)
#' ### Definition of the nature of the priors
#' type.prior=c("gaussian","gaussian","gamma")
#' ### Definition of the prior hyperparameters
#' opt.prior=list(c(11,3),c(6,2),c(2,0.1))
#' ### Definition of the estimation option
#' opt.estim=list(Ngibbs=400,Nmh=1000,thetaInit=c(11,6,0.1),k=rep(5e-4,3),sig=diag(3))
#' ### Definition of the emulation options
#' opt.emul=list(p=1,n.emul=150,PCA=FALSE)
#'
#' test <- estim(code,X,Yr,Yexp,model="model1",type.prior,log=TRUE,opt.prior,opt.estim)
#' test2 <- estim(code,X,Yr,Yexp,model="model2",type.prior,log=TRUE,opt.prior,opt.estim,opt.emul)
#' test3 <- estim(code,X,Yr,Yexp,model="model3",type.prior,log=TRUE,opt.prior,opt.estim)
#' test4 <- estim(code,X,Yr,Yexp,model="model4",type.prior,log=TRUE,opt.prior,opt.estim,opt.emul)
#' test$plot()
#' test2$plot()
#' test3$plot()
#' test4$plot()
#'
#'
#' #### with two X
#' X <- cbind(seq(0,1,length.out=200),seq(0,1,length.out=200))
#' code <- function(X,theta)
#' {
#'   return((theta[2]*X[,1]-2)^2*sin(theta[1]*X[,2]-4))
#' }
#' ### Simulated data
#' Yr   <- code(X,c(11,6))
#' Yexp <- Yr+rnorm(200,0,0.1)
#' ### Definition of the nature of the priors
#' type.prior=c("gaussian","gaussian","gamma")
#' ### Definition of the prior hyperparameters
#' opt.prior=list(c(11,3),c(6,2),c(2,0.1))
#' ### Definition of the estimation option
#' opt.estim=list(Ngibbs=400,Nmh=1000,thetaInit=c(11,6,0.1),k=rep(5e-4,3),sig=diag(3))
#'
#' test <- estim(code,X,Yr,Yexp,model="model1",type.prior,log=TRUE,opt.prior,opt.estim)
#' test2 <- estim(code,X,Yr,Yexp,model="model2",type.prior,log=TRUE,opt.prior,opt.estim,opt.emul)
#' test3 <- estim(code,X,Yr,Yexp,model="model3",type.prior,log=TRUE,opt.prior,opt.estim)
#' test4 <- estim(code,X,Yr,Yexp,model="model4",type.prior,log=TRUE,opt.prior,opt.estim,opt.emul)
#' test$plot(depend.X=FALSE)
#' test2$plot(depend.X=FALSE)
#' test3$plot(depend.X=FALSE)
#' test4$plot(depend.X=FALSE)
#'
#'
#' @export
estim <-function(code,X,Yr,Yexp,model="model1",type.prior,log=TRUE
                  ,opt.prior,opt.estim,opt.emul=list(p=1,n.emul=100,PCA=TRUE),type.valid=NULL,opt.valid=NULL)
{
  res <- estim.class$new(code,X,Yr,Yexp,model,type.prior,log,opt.emul,opt.prior,opt.estim,type.valid,opt.valid)
  return(res)
}


#' Generates \code{\link{predict.class}} objects
#'
#' \code{predict} is a function that allows us to generate a class in which the estimation is
#' done
#'
#' The realized estimation is realized similarly as it is defined in [1]
#'
#' @useDynLib calibrationCode
#'
#' @param est a \code{\link{estim.class}} object
#' @param x.new newdata for the prediction
#' @return lala for the moment
#' @author M. Carmassi
#' @seealso \code{\link{model.class}}, \code{\link{prior.class}}
#' @examples
#' ####### For the first model
#' ### The data set
#' X <- seq(0,1,length.out=100)
#' ### The code to calibrate
#' code <- function(X,theta)
#' {
#'   return((6*X-2)^2*sin(theta*X-4))
#' }
#' ### Simulated data
#' Yr   <- code(X,11)
#' Yexp <- Yr+rnorm(100,0,0.1)
#' ### Definition of the nature of the priors
#' type.prior=c("gaussian","gamma")
#' ### Definition of the prior hyperparameters
#' opt.prior=list(c(11,3),c(2,0.1))
#' ### Definition of the estimation option
#' opt.estim=list(Ngibbs=400,Nmh=1000,thetaInit=c(11,0.1),k=c(5e-4,5e-4),sig=diag(2))
#'
#' test <- estim(code,X,Yr,Yexp,model="model1",type.prior,log=TRUE,opt.prior,opt.estim)
#' test$plot()
#'
#' @export
prediction <-function(est,x.new)
{
  res <- predict.class$new(est,x.new)
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

