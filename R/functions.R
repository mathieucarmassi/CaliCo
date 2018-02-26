#' Generates \code{\link{model.class}} objects.
#'
#' \code{model} is a function that generates a calibration model and the associated likelihood.
#'
#'
#' @details The different statistical models are: \itemize{\item{Model1:
#' \deqn{for i in [1,...,n]  Yexp_i=f(x_i,\Theta)+\epsilon(x_i)}}
#' \item{Model2:
#' \deqn{for i in [1,...,n]  Yexp_i=F(x_i,\Theta)+\epsilon(x_i)}}
#' \item{Model3:
#' \deqn{for i in [1,...,n]  Yexp_i=f(x_i,\Theta)+\delta(x_i)+\epsilon(x_i)}}
#' \item{Model4:
#' \deqn{for i in [1,...,n]  Yexp_i=F(x_i,\Theta)+\delta(x_i)+\epsilon(x_i)}}
#' }
#' where \eqn{for i in [1,\dots,n] \epsilon(x_i)~N(0,\sigma^2)}, \eqn{F(.,.)~PG(m_1(.,.),c_1{(.,.),(.,.)})}
#'  and \eqn{\delta(.)~PG(m_2(.),c_2(.,.))}.
#' There is four kind of models in calibration. They are properly defined in [1].
#'
#'
#' @param code the computational code (function of X and theta)
#' @param X the matrix of the forced variables
#' @param Yexp the vector of the experiments
#' @param model string of the model chosen ("model1","model2","model3","model4")
#' by default "model1" is choosen. See details for precisions.
#' @param opt.emul is an option list containing characteristics about emulation option. \itemize{
#' \item{\strong{p}}{ the number of parameter in the model (defaul value 1)}
#' \item{\strong{n.emul}}{ the number of points for contituing the Design Of Experiments (DOE) (default value 100)}
#' \item{\strong{type}}{ type of the chosen kernel (value by default "matern5_2") from \code{\link{km}} function}
#' \item{\strong{binf}{ the lower bound of the parameter vector (default value 0)}}
#' \item{\strong{bsup}{ the upper bound of the parameter vector (default value 1)}}
#' \item{\strong{DOE}{ design of experiments for the surrogate (default value NULL). If NULL the DOE is automatically
#' generated as a maximin LHS}}
#' }
#' @param opt.disc is an option list containing characteristics on the discrepancy \itemize{
#' \item{\strong{kernel.type}{ see \code{\link{kernel.fun}} for further details}}
#' }
#' @return \code{model} returns a \code{model.class} object. This class contains two main methods:
#' \itemize{
#' \item{$plot(\eqn{\Theta},\eqn{\sigma^2}, points=FALSE)}{ this metod generates the plot for a new
#' \eqn{\Theta}, \eqn{\sigma^2} and a new set of data. The option points allows to vizualize the points from
#' the Design Of Experiments (DOE) used for establishing the surrogate.}
#' \item{$print()}{ this method presents the main information about the model.}
#' }
#' @author M. Carmassi
#' @seealso \code{\link{prior}},\code{\link{calibrate}},\code{\link{prediction}}, \code{\link{kernel.fun}}
#' @examples
#' ###### The code to calibrate
#' X <- cbind(seq(0,1,length.out=100),seq(0,1,length.out=100))
#' code <- function(X,theta)
#' {
#'   return((6*X[,1]*theta[2]-2)^2*theta[1]*sin(theta[3]*X[,2]-4))
#' }
#' Yexp <- code(X,c(1,1,11))+rnorm(100,0,0.1)
#'
#' ###### For the first model
#' ### Generate the model
#' model1 <- model(code,X,Yexp,"model1")
#' ### Plot the results with the first column of X
#' model1$plot(c(1,1,11),0.1,select.X=X[,1])
#' ### Summury of the foo class generated
#' model1$print()
#' ### Acces of the fun in the model
#' model1$fun(c(1,1,11),0.1)
#' ### Get acces to the likelihood
#' model1$likelihood(c(1,1,11),0.1)
#'
#' ###### For the second model
#' ### Generate the model with setup for the Gaussian Process
#' binf <- c(0.9,0.9,10.5)
#' bsup <- c(1.1,1.1,11.5)
#' opt.emul <- list(p=3,n.emul=100,type="matern5_2",binf=binf,bsup=bsup,DOE=NULL)
#' model2 <- model(code,X,Yexp,"model2",opt.emul)
#' ### Plot the model
#' model2$plot(c(1,1,11),0.1,select.X=X[,1])
#'
#' ### Use your own design of experiments
#' DOE <- DiceDesign::lhsDesign(100,5)$design
#' DOE[,3:5] <- unscale(DOE[,3:5],binf,bsup)
#' opt.emul <- list(p=3,n.emul=100,type="matern5_2",binf=c(0.9,0.9,10.5),bsup=c(1.1,1.1,11.5),DOE=DOE)
#' model2 <- model(code,X,Yexp,"model2",opt.emul)
#' plot(model2, theta=c(1,1,11),sig=0.1,points=FALSE,select.X=X[,1])
#' print(model2)
#'
#'
#' ###### For the third model
#' model3 <- model(code,X,Yexp,"model3",opt.disc=list(kernel.type="matern5_2"))
#' model3$plot(c(1,1,11),c(2,0.5),0.1,select.X=X[,1])
#' model3$likelihood(c(1,1,11),c(2,0.5),0.1)
#' model3$print()
#'
#' ###### For the fourth model
#' ### Desactivation of the input DOE
#' opt.emul=list(p=3,n.emul=100,type="matern5_2",binf=binf,bsup=bsup,DOE=NULL)
#' model4 <- model(code,X,Yexp,"model4",opt.emul)
#' model4$plot(c(1,1,11),c(2,0.5),0.1,select.X=X[,1])
#'
#' ### Use your own design of experiments
#' DOE <- DiceDesign::lhsDesign(100,5)$design
#' DOE[,3:5] <- unscale(DOE[,3:5],binf,bsup)
#' opt.emul <- list(p=3,n.emul=100,type="matern5_2",binf=binf,bsup=bsup,DOE=DOE)
#' model4 <- model(code,X,Yexp,"model4",opt.emul=list(p=3,n.emul=100,type="matern5_2",binf=binf,bsup=bsup,DOE=DOE))
#' model4$plot(c(1,1,11),c(2,0.5),0.1,points=FALSE,select.X=X[,1])
#' model4$print()
#' model4$likelihood(c(1,1,11),c(2,0.5),0.1)
#'
#' @export
model <- function(code,X,Yexp,model="model1",opt.emul=list(p=1,n.emul=100,type="matern5_2",
                                                           binf=0,bsup=1,DOE=NULL),opt.disc=list(kernel.type=NULL))
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
           obj = model3.class$new(code,X,Yexp,model,opt.disc)
           return(obj)
         },
         model4={
           obj = model4.class$new(code,X,Yexp,model,opt.emul,opt.disc)
           return(obj)
         }
  )
}


#' Generates \code{\link{prior.class}} objects.
#'
#' \code{prior} is a function that generates a \code{\link{prior.class}} containing information about one or
#' several priors. When several priors are selected, the function \code{prior} return a list of \code{\link{prior.class}}.
#'
#' @details The densities implemented are defined as follow
#' \itemize{
#' \item{The Gaussian density:
#' \deqn{f(x)=1/(\sigma*\sqrt(2\pi))exp{-1/2*((x-\mu)/\sigma)^2}}
#' where \strong{\eqn{\mu}} and \strong{\eqn{\sigma}} (the mean and the standard deviation)
#' are the two hyperparameters. The vector \eqn{c(\mu,\sigma)} is the one looked for in opt.prior.}
#' \item{The Gamma density:
#' \deqn{f(x)=1/(k^a*\Gamma(a))*x^(a-1)*exp(-(x/k))}
#' where \strong{\eqn{a}} and \strong{\eqn{k}} (the shape and the scale)
#' are the two hyperparameters. The vector \eqn{c(a,k)} is the one looked for in opt.prior.}
#' \item{The Uniform density:
#' \deqn{f(x)=1/(b-a)}
#' where \strong{\eqn{a}} and \strong{\eqn{b}} (the upper and the lower bound)
#' are the two hyperparameters. The vector \eqn{c(a,b)} is the one looked for in opt.prior.}
#' }
#'
#' @import ggplot2
#' @import gridExtra
#' @param  type.prior the vector of the prior types selected. For example type.prior=c("gaussian","beta")
#' @param opt.prior list of the hyperparameters relatives to the prior selected. If the first prior selected is
#' Gaussian, the hyperparameters would be the mean and the standard deviation. See Details for precisions.
#' @param log (default=TRUE) if the log value is wanted or not.
#' @return \code{prior} returns a \code{\link{prior.class}} object. Two main methods are available:
#' \itemize{\item{$plot()}{ display the probability density of the prior}
#' \item{$print()}{ return the main information concerning the prior distribution}}
#' @author M. Carmassi
#' @seealso \code{\link{prior}},\code{\link{calibrate}},\code{\link{prediction}}, \code{\link{kernel.fun}}
#' @examples
#' #### Only one prior is wanted
#' ## For a Gaussian Prior
#' gaussian <- prior(type.prior="gaussian",opt.prior=list(c(0.5,0.001)))
#' gaussian$plot()
#'
#' ## For a Uniform Prior
#' unif <- prior(type.prior="unif",opt.prior=list(c(0,1)))
#' unif$plot()
#'
#' ## For a Gamma Prior
#' gamma <- prior(type.prior="gamma",opt.prior=list(c(5,1)))
#' gamma$plot()
#'
#' #### For several priors
#' priors <- prior(type.prior=c("gaussian","gamma"),opt.prior=list(c(0.5,0.001),c(5,1)))
#' grid.arrange(priors$Prior1$plot(),priors$Prior2$plot(),nrow=2)
#'
#' @export
prior <- function(type.prior,opt.prior,log=TRUE)
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


#' Generates \code{\link{calibrate.class}} objects
#'
#' \code{calibration} is a function that allows us to generate a \code{\link{calibrate.class}} class in which the estimation is
#' done from a \code{\link{model.class}} and a \code{\link{prior.class}} objects.
#'
#' @useDynLib CaliCo
#' @importFrom Rcpp evalCpp
#'
#' @param md a \code{\link{model.class}} object
#' @param pr a \code{\link{prior.class}} object
#' @param opt.estim estimation options \itemize{\item{Ngibbs}{Number of iteration of the algorithm Metropolis within Gibbs}
#' \item{Nmh}{ Number of iteration of the Metropolis Hastings algorithm}
#' \item{thetaInit}{ Initial point}
#' \item{k}{ Tuning parameter for the covariance matrix sig}
#' \item{sig}{ Covariance matrix for the proposition distribution (\eqn{k*sig})}}
#' @param onlyCV if TRUE run the cross validation only (default onlyCV=FALSE)
#' @param opt.valid list of cross validation options (default opt.valid=FALSE)\itemize{
#' \item{nCV}{ Number of iterations for the cross validation}
#' \item{type.valid}{ Type of cross validation selected. "loo" (leave one out) is the only method emplemented so far.}
#' }
#' @return \code{calibrate} returns a \code{\link{calibrate.class}} object. Two main methods are available:
#' \itemize{\item{$plot()}{ display the probability density of the prior with different options:}
#' \itemize{
#' \item {graph}{ The vector of the graph wanted. By default all the graph are displayed and graph=c("acf","chains","densities","output").
#' "acf" displays the correlation graph of the MCMC chains, "chains" plot the chains, "densities" shows the comparison of the
#' densities a priori and a posteriori, and "output" displays the output of the code with the calibrated one and its credibility
#' interval (if CI=TRUE).}
#' \item {separated}{ Allows to separate each graphs by displying one by one all the graphs. By default separated=FALSE}
#' \item {CI}{ Allows to add the posterior credibility interval to the output plot. By default CI=TRUE}
#' \item {select.X}{ When the number of X is >1, this option has to be activated to display the output plot. select.X
#' allows to choose one X for the x scale in the output plot}}
#' \item{$print()}{ return the main information concerning the estim.class object}}
#' @author M. Carmassi
#' @seealso \code{\link{prior}},\code{\link{calibrate}},\code{\link{prediction}}, \code{\link{kernel.fun}}
#' @examples
#' ###################### The code to calibrate
#' X <- cbind(seq(0,1,length.out=10),seq(0,1,length.out=10))
#' code <- function(X,theta)
#' {
#'   return((6*X[,1]*theta[2]-2)^2*theta[1]*sin(theta[3]*X[,2]-4))
#' }
#' Yexp <- code(X,c(1,1,11))+rnorm(10,0,0.1)
#'
#' ############### For the first model
#' ###### Definition of the model
#' md <- model(code,X,Yexp,"model1")
#' ###### Definition of the prior densities
#' pr <- prior(type.prior=c("gaussian","gaussian","gaussian","gamma"),opt.prior=
#' list(c(1,0.01),c(1,0.01),c(11,3),c(2,0.1)))
#' ###### Definition of the calibration options
#' opt.estim=list(Ngibbs=2000,Nmh=6000,thetaInit=c(1,1,11,0.1),k=c(6e-3,1e-3,1e-5,1e-3),
#' sig=diag(4),Nchains=1,burnIn=3000)
#' ###### Run the calibration
#' mdfit <- calibrate(md,pr,opt.estim)
#' ####### The plot generated is a list of ggplot
#' p <- plot(mdfit,select.X=X[,1])
#' print(mdfit)
#' ####### Run regular calibration plus cross validation (the plot function is available)
#' opt.valid <- list(type.valid='loo',nCV=10)
#' mdfitCV <- calibrate(md,pr,opt.estim,opt.valid)
#' print(mdfitCV)
#' ####### Run cross validataion only (the plot function is disabled)
#' mdfitCV2 <- calibrate(md,pr,opt.estim,opt.valid,onlyCV=TRUE)
#' print(mdfitCV2)
#'
#'
#' ############### For the second model
#' ###### Definition of the model
#' ### The lower and upper bound vector for the parameter
#' binf <- c(0.9,0.9,10.5)
#' bsup <- c(1.1,1.1,11.5)
#' opt.emul <- list(p=3,n.emul=50,type="matern5_2",binf=binf,bsup=bsup,DOE=NULL)
#' md2 <- model(code,X,Yexp,"model2",opt.emul)
#' ###### Definition of the prior densities
#' pr <- prior(type.prior=c("gaussian","gaussian","gaussian","gamma"),opt.prior=
#' list(c(1,0.01),c(1,0.01),c(11,3),c(2,0.1)))
#' ###### Definition of the calibration options
#' opt.estim=list(Ngibbs=2000,Nmh=6000,thetaInit=c(1,1,11,0.1),k=c(6e-3,1e-3,1e-5,1e-3),
#' sig=diag(4),Nchains=1,burnIn=3000)
#' ###### Run the calibration
#' mdfit2 <- calibrate(md2,pr,opt.estim)
#' ####### The plot generated is a list of ggplot
#' p <- plot(mdfit2,select.X=X[,1])
#' print(mdfit2)
#' ####### Run cross validataion only (the plot function is disabled)
#' opt.valid <- list(type.valid='loo',nCV=10)
#' mdfitCV2 <- calibrate(md2,pr,opt.estim,opt.valid,onlyCV=TRUE)
#' print(mdfitCV2)
#'
#'
#' ### For the third model
#' md3 <- model(code,X,Yexp,"model3",opt.disc=list(kernel.type="gauss"))
#'
#' ### For the fourth model
#' md4 <- model(code,X,Yexp,"model4",opt.emul,opt.disc=list(kernel.type="matern5_2"))
#'
#' ### Definition of the priors
#' pr1 <- prior(type.prior=c("gaussian","gaussian","gaussian","gamma"),opt.prior=
#' list(c(1,0.01),c(1,0.01),c(11,3),c(2,0.1)))
#'
#' pr2 <- prior(type.prior=c("gaussian","gaussian","gaussian","gaussian","gamma","gamma"),opt.prior=
#' list(c(1,0.01),c(1,0.01),c(11,3),c(2,0.1),c(2,0.1),c(2,0.1)))
#'
#' ### Calibration with estimation options
#' opt.estim1=list(Ngibbs=2000,Nmh=6000,thetaInit=c(1,1,11,0.1),k=c(6e-3,1e-3,1e-5,1e-3),sig=diag(4),Nchains=1,burnIn=3000)
#' opt.estim2=list(Ngibbs=2000,Nmh=6000,thetaInit=c(1,1,11,2,0.1,0.1),k=rep(5e-3,6),sig=diag(6),Nchains=1,burnIn=3000)
#'
#' ### Calibration model1
#' modelfit <- calibrate(md1,pr1,opt.estim1)
#' p <- modelfit$plot(select.X=X[,1])
#' modelfit$print()
#' # Calibration plus cross validation
#' opt.valid <- list(type.valid='loo',nCV=4)
#' modelfitCV <- calibrate(md1,pr1,opt.estim1,opt.valid)
#' p <- modelfitCV$plot(select.X=X[,1])
#' modelfitCV$print()
#' # Cross validataion only
#' modelfitCV2 <- calibrate(md1,pr1,opt.estim1,opt.valid,onlyCV=TRUE)
#' modelfitCV2$print()
#'
#' ### Calibration model2
#' modelfit2 <- calibrate(md2,pr1,opt.estim1)
#' p <- modelfit2$plot(select.X=X[,1])
#' opt.valid <- list(type.valid='loo',nCV=4)
#' # Cross validataion only
#' modelfitCV2 <- calibrate(md2,pr1,opt.estim1,opt.valid,onlyCV=TRUE)
#' p <- modelfit2$plot(select.X=X[,1])
#'
#' ### Calibration model3
#' modelfit3 <- calibrate(md3,pr2,opt.estim2)
#' opt.valid <- list(type.valid='loo',nCV=4)
#' # Cross validataion only
#' modelfit3CV <- calibrate(md3,pr2,opt.estim2,opt.valid,onlyCV=TRUE)
#' p <- modelfit3$plot(select.X=X[,1])
#'
#' ### Calibration model4
#' modelfit4 <- calibrate(md4,pr2,opt.estim2)
#' # Cross validataion only
#' modelfit4CV <- calibrate(md4,pr2,opt.estim2,opt.valid,onlyCV=TRUE)
#' p <- modelfit4$plot(select.X=X[,1])
#'
#' @export
calibrate <-function(md,pr,opt.estim,opt.valid=NULL,onlyCV=FALSE)
{
  res <- calibrate.class$new(md,pr,opt.estim,opt.valid,onlyCV)
  return(res)
}


#' Generates \code{\link{predict.class}} objects
#'
#' \code{predict} is a function that allows us to generate a class in which the estimation is
#' done
#'
#' The realized estimation is realized similarly as it is defined in [1]
#'
#' @useDynLib CaliCo
#'
#' @param modelfit a \code{\link{estim.class}} object
#' @param newdata newdata for the prediction
#' @return return a \code{\link{predict.class}} object with two main methods
#' @author M. Carmassi
#' @seealso \code{\link{model.class}}, \code{\link{prior.class}}, \code{\link{esim.class}}
#' @examples
#' ### The code to calibrate
#' X <- cbind(seq(0,1,length.out=10),seq(0,1,length.out=10))
#' code <- function(X,theta)
#' {
#'   return((6*X[,1]*theta[2]-2)^2*theta[1]*sin(theta[3]*X[,2]-4))
#' }
#' Yexp <- code(X,c(1,1,11))+rnorm(10,0,0.1)
#'
#'
#' # Definition of the different models
#' md1 <- model(code,X,Yexp,"model1")
#'
#' ### For the second model
#' binf <- c(0.9,0.9,10.5)
#' bsup <- c(1.1,1.1,11.5)
#' DOE <- DiceDesign::lhsDesign(10,5)$design
#' DOE[,3:5] <- unscale(DOE[,3:5],binf,bsup)
#' DOE[,1:2] <- unscale(DOE[,1:2],c(0,0),c(1.5,1.5))
#' opt.emul <- list(p=3,n.emul=10,type="matern5_2",binf=binf,bsup=bsup,DOE=DOE)
#' md2 <- model(code,X,Yexp,"model2",opt.emul)
#'
#' ### For the third model
#' md3 <- model(code,X,Yexp,"model3",opt.disc=list(kernel.type="gauss"))
#'
#' ### For the fourth model
#' md4 <- model(code,X,Yexp,"model4",opt.emul,opt.disc=list(kernel.type="matern5_2"))
#'
#' ### Definition of the priors
#' pr1 <- prior(type.prior=c("gaussian","gaussian","gaussian","gamma"),opt.prior=
#' list(c(1,0.01),c(1,0.01),c(11,3),c(2,0.1)))
#'
#' pr2 <- prior(type.prior=c("gaussian","gaussian","gaussian","gaussian","gamma","gamma"),opt.prior=
#' list(c(1,0.01),c(1,0.01),c(11,3),c(2,0.1),c(2,0.1),c(2,0.1)))
#'
#' ### Calibration with estimation options
#' opt.estim1=list(Ngibbs=400,Nmh=600,thetaInit=c(1,1,11,0.1),k=rep(5e-4,4),sig=diag(4),Nchains=1,burnIn=300)
#' opt.estim2=list(Ngibbs=400,Nmh=600,thetaInit=c(1,1,11,2,0.1,0.1),k=rep(5e-3,6),sig=diag(6),Nchains=1,burnIn=300)
#'
#' ### Calibration model1
#' modelfit <- calibrate(md1,pr1,opt.estim1)
#'
#' ###
#' x.new <- cbind(seq(1,1.5,length.out=10),seq(1,1.5,length.out=10))
#' emul <- prediction(modelfit,x.new)
#' emul$plot(select.X=x.new[,1])
#'
#'
#' modelfit2 <- calibrate(md2,pr1,opt.estim1)
#' ###Prediction
#' x.new <- cbind(seq(1,1.5,length.out=10),seq(1,1.5,length.out=10))
#' emul <- prediction(modelfit2,x.new)
#' emul$plot(select.X=x.new[,1])
#'
#' ### For the third model
#' modelfit3 <- calibrate(md3,pr2,opt.estim2)
#' ### Prediction
#' x.new <- cbind(seq(1,1.5,length.out=10),seq(1,1.5,length.out=10))
#' emul <- prediction(modelfit3,x.new)
#' emul$plot(select.X=x.new[,1])
#'
#'
#' @export
prediction <-function(modelfit,x.new)
{
  res <- prediction.class$new(modelfit,x.new)
  return(res)
}


#' Generates covariances matrices thanks to \code{\link{Kernel.class}}
#'
#' \code{Kernel.fun} is a function that allows us to generate covariances matrices from data
#'
#' @param X data
#' @param var the variance for the covariance function
#' @param psi the parameter vector
#' @param kernel.type the choice of the form of the kernel (with d chosen as an euclidian distance) \itemize{
#' \item{gauss}{ \deqn{\sigma^2 exp{-1/2(d/\psi)^2}}}
#' \item{exp}{ \deqn{\sigma^2 exp{-1/2 d/\psi}}}
#' \item{matern3_2}{ \deqn{\sigma^2(1+\sqrt{3}d^2/\psi) exp{-\sqrt{3}d^2/\psi}}}
#' \item{matern5_2}{ \deqn{\sigma^2(1+\sqrt{5}d^2/\psi+5d^2/(3\psi^2))exp{-\sqrt{5}d^2/\psi}}}}
#' @return \code{Kernel.fun} returns a covariance matrix
#' @author M. Carmassi
#' @seealso \code{\link{model.class}}, \code{\link{prior.class}}, \code{\link{esim.class}}
#' @examples
#' X <- cbind(seq(0,10,length.out=10),seq(8,20,length.out=10))
#' var <- 2
#' psi <- 0.1
#' Cov <- kernel.fun(X,var,psi,kernel.type="matern5_2")
#'
#' @export
kernel.fun <- function(X,var,psi,kernel.type="gauss")
{
  library(R6)
  if (is.null(kernel.type)){kernel.type <- "gauss"}
  switch(kernel.type,
         gauss={
           obj = gauss.class$new(X,var,psi,kernel.type)
           return(obj$Cov)
         },
         exp={
           obj = exp.class$new(X,var,psi,kernel.type)
           return(obj$Cov)
         },
         matern3_2={
           obj = matern3_2.class$new(X,var,psi,kernel.type)
           return(obj$Cov)
         },
         matern5_2={
           obj = matern5_2.class$new(X,var,psi,kernel.type)
           return(obj$Cov)
         }
  )
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


#' Function that deals with negative eigen values in a matrix not positive definite
#'
#' @param M the matrix or the vector
#' @return the new matrix M
#' @export
DefPos <- function(X)
{
  p <- eigen(X)$vectors
  e <- eigen(X)$values
  if (all(e>0)){} else
  {
    e[which(e<0)] <- 1e-4
  }
  d <- diag(e)
  return(t(p)%*%d%*%p)
}
