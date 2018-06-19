#' A Reference Class to generates differents model objects
#'
#'
#' @description See the function \code{\link{model}} which produces an instance of this class
#' This class comes with a set of methods, some of them being useful for the user:
#' See the documentation for \code{\link{model}}... Other methods
#' should not be called as they are designed to be used during the calibration process.
#'
#'
#' Fields should not be changed or manipulated by the user as they are updated internally
#' during the estimation process.
#' @field code a function which takes in entry X and theta
#' @field X the matrix of the forced variables
#' @field Yexp the experimental output
#' @field n the number of experiments
#' @field d the number of forced variables
#' @field binf the lower bound of the parameters for the DOE
#' @field bsup the upper bound of the parameters for the DOE
#' @field opt.gp a list of parameter for the surrogate (default NULL) \itemize{
#' \item{\strong{type}}{ type of the chosen kernel (value by default "matern5_2") from \code{\link{km}} function}
#' \item{\strong{DOE}}{ design of experiments for the surrogate (default value NULL)}}
#' @field opt.emul a list of parameter to establish the DOE (default NULL) \itemize{
#' \item{\strong{p}}{ the number of parameter in the model (default value 1)}
#' \item{\strong{n.emul}}{ the number of points for constituting the Design Of Experiments (DOE) (default value 100)}
#' \item{\strong{binf}}{ the lower bound of the parameter vector (default value 0)}
#' \item{\strong{bsup}}{ the upper bound of the parameter vector (default value 1)}}
#' @field opt.sim a list of parameter containing output of the code and corresponding DOE \itemize{
#' \item{\strong{Ysim}}{ Output of the code}
#' \item{\strong{DOEsim}}{ DOE corresponding to the output of the code}}
#' @field model the model choice (see \code{\link{model}} for more specification).
#' @field opt.disc a list of parameter for the discrepancy \itemize{
#' \item{\strong{kernel.type}}{ the kernel choosen for the Gaussian process}}
#'
#' @export
model.class <- R6Class(classname = "model.class",
                 public = list(
                   code      = NULL, ## code variable
                   X         = NULL, ## forced variables
                   Yexp      = NULL, ## experiments
                   n         = NULL, ## length of the experiements
                   d         = NULL, ## number of forced variables
                   model     = NULL, ## statistical model elected
                   theta     = NULL, ## parameter vector
                   var       = NULL, ## variance of the measurement error
                   thetaD    = NULL, ## discrepancy parameter vector
                   n.cores   = NULL, ## number of computer cores
                   initialize = function(code=NA,X=NA,Yexp=NA,model=NA)
                   {
                     self$code    <- code
                     self$X       <- as.matrix(X)
                     self$Yexp    <- Yexp
                     self$n       <- length(Yexp)
                     self$model   <- model
                     self$n.cores <- detectCores()
                     if (is.matrix(X)) {self$d <- ncol(X)} else{self$d <-1}
                     private$checkModels()
                     private$checkCode()
                     private$checkOptions()
                   },
                   active = list(
                     theta  = function(value) {return(value)},
                     thetaD = function(value) {return(value)},
                     var    = function(value) {return(value)}
                   ),
                   gglegend =function()
                   {
                     return(theme(legend.title = element_blank(),
                                  legend.position = c(0.2,0.8),
                                  legend.background = element_rect(
                                    linetype="solid", colour ="grey")))
                   }
                 ))

model.class$set("private","checkModels",
        function()
          ### Check if the chosen model is in the possible selection
        {
          if (!self$model %in% c("model1","model2","model3","model4"))
          {
            stop('Please elect a correct model')
          }
        })


model.class$set("private","checkOptions",
                function()
                  ### Check if there is no missing in the options
                  {
                  test <- function(N,options)
                  {
                    N2 <- names(options)
                    for (i in 1:length(N))
                    {
                      if(names(options)[i] != N[i])
                      {
                        stop(paste(N[i],"value is missing, please enter a correct value",sep=" "))
                      }
                    }
                  }
                  if (self$model %in% c("model2","model4"))
                  {
                    if (!is.null(self$opt.emul)) test(c("p","n.emul","binf","bsup"),self$opt.emul)
                    test(c("type","DOE"),self$opt.gp)
                    if (!is.null(self$opt.sim)) test(c("Ysim","DOEsim"),self$opt.sim)
                  } else
                  {
                    if (self$model %in% c("model3","model4")){test(c("kernel.type"),self$opt.disc)}
                  }
                })

model.class$set("private","checkCode",
                function()
                  ### Check if the code is valid
                {
                  if (is.null(self$code) & self$model %in% c("model1","model3"))
                  {
                    stop("The code cannot be NULL if you chose model1 or model3")
                  }
                })


################################## Model 1 definition #######################################

model1.class <- R6Class(classname = "model1.class",
                        inherit = model.class,
                        public=list(
                        m.exp = NULL, ## Mean for the likelihood
                        V.exp = NULL, ## Variance for the likelihood
                        l     = NULL, ## Length of the new data for prediction
                        initialize=function(code=NA, X=NA, Yexp=NA, model=NA)
                        {
                          ### Initialize from model.class
                          super$initialize(code, X, Yexp, model)
                        },
                        model.fun = function(theta,var,X=self$X,CI="err")
                        {
                          ### Function that generates the output of the model. If CI=TRUE, it computes the credibility
                          ### intervals of the white Gaussian noise
                          y <- matrix(nr=100,nc=self$n)
                          for(i in 1:100){y[i,] <- self$code(X,theta)+rnorm(self$n,0,sqrt(var))}
                          qq <- apply(y,2,quantile,c(0.05,0.5,0.95))
                          if (is.null(CI))
                          {
                            df <- data.frame(y=qq[2,],type="model output")
                          } else if (CI=="err")
                          {
                            df <- data.frame(y=qq[2,],type="model output",q05=qq[1,],q95=qq[3,],fill="CI 90% noise")
                          } else
                          {
                            warning("The argument for the credibility interval is not valid and no credibility interval
                                    will be displayed")
                            df <- data.frame(y=qq[2,],type="model output")
                          }
                          return(df)
                        }
                        # prediction.fun = function(theta,var,x.new)
                        # {
                        #   ### Prediction function is the function to use when applying on a new data set
                        #   if (is.matrix(x.new)){l <- nrow(x.new)} else{l <- length(x.new)}
                        #   y  <- self$code(x.new,theta)+rnorm(l,0,sqrt(var))
                        #   df <- data.frame(pr=y,type="predicted")
                        #   return(df)
                        #}
                        )
)

## likelihood function
model1.class$set("public","likelihood",
                function(theta,var)
                {
                  ### Log-Likelihood
                  self$m.exp <- self$code(self$X,as.vector(theta))
                  self$V.exp <- var*diag(self$n)
                  return(-self$n/2*log(2*pi)-1/2*log(det(self$V.exp))
                         -0.5*t(self$Yexp-self$m.exp)%*%solve(self$V.exp)%*%(self$Yexp-self$m.exp))
                })

## print function
model1.class$set("public","print",
                 function(...)
                 {
                   cat("Call:\n")
                   print(self$model)
                   cat("\n")
                   cat("With the function:\n")
                   print(self$code)
                   cat("\n")
                   cat("No surrogate is selected")
                   cat("\n")
                   cat("No discrepancy is added")
                 })

## plot function
model1.class$set("public","plot",
                 function(x,CI="err",...)
                 {
                   ### Plot generates a ggplot object
                   if (is.matrix(x)){stop("please enter a correct x to plot your model")}
                   if (length(x)!= self$n){stop(paste("please enter a correct vector x of size",
                                                      self$n,sep=" "))}
                   df <- data.frame(y=self$Yexp,type="exp")
                   if (!is.null(self$theta) & !is.null(self$var)) {
                     df2 <- self$model.fun(self$theta,self$var,self$X,CI)
                     df2 <- cbind(df2,x=x)
                     if (is.null(CI))
                     {
                       df  <- cbind(df,x=x)
                       df  <- rbind(df,df2)
                       p <- ggplot(df) + geom_line(mapping = aes(x=x,y=y,color=type))+
                         xlab("")+ylab("")+theme_light()+self$gglegend()+
                         scale_color_manual(values=c("red", "#000000"))
                     } else if (CI=="err")
                     {
                       df  <- cbind(df,q05=df2$q05,q95=df2$q95,fill="CI 90% noise",x=x)
                       df  <- rbind(df,df2)
                       p <- ggplot(df) + geom_line(mapping = aes(x=x,y=y,color=type))+
                         xlab("")+ylab("")+theme_light()+self$gglegend()+
                         geom_ribbon(mapping = aes(x=x,ymin=q05,ymax=q95,fill=fill),alpha=0.4)+
                         scale_fill_manual(values = "skyblue3")+
                         scale_color_manual(values=c("red", "#000000"))
                     } else {
                       df  <- cbind(df,x=x)
                       df  <- rbind(df,df2)
                       p <- ggplot(df) + geom_line(mapping = aes(x=x,y=y,color=type))+
                         xlab("")+ylab("")+theme_light()+self$gglegend()
                       scale_color_manual(values=c("red", "#000000"))
                     }
                   } else
                   {
                     p <- ggplot(df) + geom_line(mapping = aes(x=x,y=y,color=type))+
                       xlab("")+ylab("")+theme_light()+self$gglegend()+
                       scale_color_manual(values=c("red"))
                   }
                   return(p)
                 })


##################################### Model 3 definition ##################################

## model main functions
model3.class <- R6Class(classname = "model3.class",
                        inherit = model1.class,
                        public=list(
                          model1.fun        = NULL, ## model function from model1
                          model1.prediction = NULL, ## prediction function from model1
                          opt.disc          = NULL, ## discrepancy options
                          disc              = NULL, ## discrepancy field
                          initialize=function(code=NA, X=NA, Yexp=NA, model=NA,opt.disc=list(kernel.type=NULL))
                          {
                            ## Check if the opt.emul option is filled if it is not a gaussian kernel is picked
                            if (is.null(opt.disc$kernel.type))
                            {
                              warning("default value is selected. The discrepancy will have a Matérn 5/2 covariance
                                      structure")
                              self$opt.disc$kernel.type="matern5_2"
                            } else
                            {
                              self$opt.disc  <- opt.disc
                            }
                            ## Initialize with the model1.class which initialize from model.class
                            super$initialize(code, X, Yexp, model)
                            ## Store the model1 functions
                            self$model1.fun        <- super$model.fun
                            self$model1.prediction <- super$prediction.fun
                          },
                          ## Parallel computation for the discrepancy
                          disc.par = function(i)
                          {
                            return(self$model1.fun(self$theta,self$var,self$X,CI=NULL)$y+
                                     self$discrepancy(self$theta,self$thetaD,self$var,self$X)$bias)
                          },
                          model.fun = function(theta,thetaD,var,X=self$X,CI="err")
                          {
                            res.model1 <- self$model1.fun(theta,var,X,CI=CI)
                            self$disc  <- self$discrepancy(theta,thetaD,var,X)
                            if (is.null(CI))
                            {
                              df <- data.frame(y=res.model1$y+self$disc$bias,
                                               type="model output")
                            } else if (CI=="err")
                            {
                              disc.temp <- do.call(rbind,mclapply(c(1:100),self$disc.par,mc.cores = self$n.cores))
                              qq <- apply(disc.temp,2,quantile,c(0.05,0.5,0.95))
                              df <- data.frame(y=qq[2,],type="model output",q05=qq[1,],q95=qq[3,],
                                               fill="CI 90% discrepancy + noise")
                            } else
                            {
                              warning("The argument for the credibility interval is not valid and no credibility interval
                                    will be displayed")
                              df <- data.frame(y=res.model1$y+qq[2,],type="model output")
                            }
                            return(df)
                          }
                          # prediction.fun = function(theta,thetaD,var,x.new)
                          # {
                          #   self$disc <- self$discrepancy(theta,thetaD,var,x.new)
                          #   foo <- self$model1.fun(theta,var,x.new)
                          #   y <- foo$y
                          #   df <- data.frame(y=self$disc$bias+y,cov=self$disc$cov,type="predicted")
                          #   return(df)
                          # }
                          )
)

## likelihood function
model3.class$set("public","likelihood",
                 function(theta,thetaD,var)
                 {
                   # Log-Likelihood
                   self$m.exp <- self$code(self$X,as.vector(theta))
                   self$V.exp <- var*diag(self$n) + self$disc$cov
                   return(-self$n/2*log(2*pi)-1/2*log(det(self$V.exp))
                          -0.5*t(self$Yexp-self$m.exp)%*%solve(self$V.exp)%*%(self$Yexp-self$m.exp))
                 })

## discrepancy function
model3.class$set("public","discrepancy",
                 ## Define method that generates a discrepancy
                 function(theta,thetaD,var,X=self$X)
                 {
                   y   <- self$model1.fun(theta,var,X)$y
                   z   <- self$Yexp - y
                   ## Compute the discrepancy covariance
                   Cov <- kernel.fun(X,thetaD[1],thetaD[2],self$opt.disc$kernel.type)
                   if (is.vector(X) & length(X)==1)
                   {} else
                   {
                     p <- eigen(Cov)$vectors
                     e <- eigen(Cov)$values
                     if (all(e>0)){} else
                     {
                       e[which(e<0)] <- 1e-4
                     }
                     d <- diag(e)
                     if (nrow(p) == 1 & ncol(p) == 1)
                     {
                       Cov <- as.numeric(p)^2*d
                     } else
                     {
                       Cov <- t(p)%*%d%*%p
                     }
                   }
                   if (is.vector(X)){long <- length(X)}else{
                     long <- nrow(X)}
                   if (long==1)
                   {
                     if (nrow(p) == 1 & ncol(p) == 1)
                     {
                       bias <- rnorm(1,0,sqrt(Cov))
                     } else
                     {
                       bias <- rnorm(n=self$n,0,sqrt(Cov))
                     }
                   } else
                   {
                     bias <- mvrnorm(100,rep(0,long),Cov)
                     dim(bias) <- c(long,100)
                     bias <- apply(bias,1,mean)
                   }
                   return(list(bias=bias,cov=Cov))
                 })

## plot function
model3.class$set("public","plot",
                 function(x,CI="all",...)
                {
                  ### Plot generates a ggplot object
                  if (is.matrix(x)){stop("please enter a correct x to plot your model")}
                  if (length(x)!= self$n){stop(paste("please enter a correct vector x of size",
                                                     self$n,sep=" "))}
                  df <- data.frame(y=self$Yexp,type="exp")
                  if (!is.null(self$theta) & !is.null(self$var)) {
                    df2 <- cbind(self$model.fun(self$theta,self$thetaD,self$var,self$X,CI),x=x)
                    if (is.null(CI))
                    {
                      df  <- cbind(df,x=x)
                      df  <- rbind(df,df2)
                      p <- ggplot(df) + geom_line(mapping = aes(x=x,y=y,color=type))+
                        xlab("")+ylab("")+theme_light()+self$gglegend()+
                        scale_color_manual(values=c("red", "#000000"))
                    } else if (CI == "err")
                    {
                      df   <- cbind(df,q05=df2$q05,q95=df2$q95,fill=df2$fill,x=x)
                      df   <- rbind(df,df2)
                      p <- ggplot(df)+geom_line(mapping = aes(x=x,y=y, color=type))+
                        xlab("")+ylab("")+theme_light()+self$gglegend()+
                        geom_ribbon(mapping = aes(x=x,ymin=q05,ymax=q95,fill="CI 90% discrepancy + noise"),
                                    alpha=0.4)+scale_fill_manual(values = adjustcolor("skyblue3"))+
                        scale_color_manual(values=c("red", "#000000"))
                    } else
                    {
                      warning("The argument for the credibility interval is not valid and no credibility interval
                                                        will be displayed")
                      df  <- cbind(df,x=x)
                      df  <- rbind(df,df2)
                      p <- ggplot(df) + geom_line(mapping = aes(x=x,y=y,color=type))+
                        xlab("")+ylab("")+theme_light()+self$gglegend()+
                        scale_color_manual(values=c("red", "#000000"))
                    }
                  } else
                  {
                    p <- ggplot(df) + geom_line(mapping = aes(x=x,y=y,color=type))+
                      xlab("")+ylab("")+theme_light()+self$gglegend()+
                      scale_color_manual(values=c("red"))
                  }
                  return(p)
                })

## print function
model3.class$set("public","print",
                 function(...)
                 {
                   if (is.null(self$theta) | is.null(self$thetaD) | is.null(self$var)){} else
                   {
                     self$disc  <- self$discrepancy(self$theta,self$thetaD,self$var,self$X)
                     bias       <- summary(self$disc$bias)
                   }
                   cat("Call:\n")
                   print(self$model)
                   cat("\n")
                   cat("With the function:\n")
                   print(self$code)
                   cat("\n")
                   cat("No surrogate is selected")
                   cat("\n\n")
                   cat("A discrepancy is added:\n")
                   cat("Summary of the bias mean:\n")
                   print(bias)
                   cat(paste("\nCovariance of the bias:",round(mean(self$disc$cov),3),"\n\n",sep=" "))
                   cat(paste("Kernel chossen: ",self$opt.disc$kernel.type,sep=""))
                 })


################################## Model 2 definition ####################################

## model 2 main functions
model2.class <- R6Class(classname = "model2.class",
                        inherit = model.class,
                        public = list(
                          opt.emul = NULL, ## DOE creation options
                          opt.sim  = NULL, ## options if the user possess the design and the output
                          opt.gp   = NULL, ## GP options
                          case     = NULL, ## case wanted by the user (depending on options)
                          doe      = NULL, ## DOE used for the surrogate
                          z        = NULL, ## output of the code for the DOE
                          GP       = NULL, ## current Gaussian process emulated
                          cov.gp   = NULL, ## covariance matrice of the gaussian process
                          yc       = NULL, ## surrogate output without noise
                          p        = NULL, ## number of parameters
                        initialize = function(code=NA, X=NA, Yexp=NA, model=NA,...)
                        {
                          if (!exists("opt.emul")){self$opt.emul <- NULL} else{self$opt.emul <- opt.emul}
                          if (!exists("opt.sim")){self$opt.sim <- NULL} else{self$opt.sim <- opt.sim}
                          if (!exists("opt.gp")){self$opt.gp <- NULL} else{self$opt.gp <- opt.gp}
                          if (!is.null(self$opt.gp$DOE)) self$opt.emul <- NULL
                          if (!is.null(self$opt.sim$DOEsim)) self$opt.emul <- NULL
                          ## Select the case wanted by the user
                          if (!is.null(self$opt.gp))
                          {
                            if (is.null(self$opt.gp$DOE) & !is.null(self$opt.emul))
                            {
                              self$case <- "1"
                              self$p <- self$opt.emul$p
                            } else if (!is.null(self$opt.gp$DOE) & is.null(self$opt.emul))
                            {self$case <- "2"} else if (!is.null(self$opt.sim)){self$case <- "3"}
                            else{stop("please enter correct options to establish the second model")}
                          }
                          ## initialize from model.class
                          super$initialize(code, X, Yexp, model)
                          if (is.null(self$code))
                          {
                            if (is.null(opt.sim))
                            {
                              print("The numerical code is desabled, please fill the opt.sim option")
                            }
                          }
                          self$GP     <- self$surrogate()
                          print("The surrogate has been set up, you can now use the function")
                        },
                        ## model function
                        model.fun = function(theta,var,X=self$X,CI="all")
                        {
                          X  <- as.matrix(X)
                          if (ncol(X) != self$d){stop("please enter a correct X")}
                          D  <- cbind(X,matrix(rep(theta,rep(nrow(X),self$p)),
                                            nr=nrow(X),nc=self$p))
                          pr <- predict(self$GP,newdata=as.data.frame(D),type="UK",
                                        cov.compute=TRUE,interval="confidence",checkNames=FALSE)
                          if (is.null(CI))
                          {
                            df <- data.frame(y=pr$mean+rnorm(nrow(X),0,sqrt(var)),
                                             type="model output")
                          } else if (CI == "all" | CI == "err")
                          {
                            nugget <- mvrnorm(n=100,pr$mean,diag(var,100))
                            qq <- apply(nugget,2,quantile,c(0.05,0.5,0.95))
                            if (CI == "all")
                            {
                              df  <- data.frame(y=qq[2,],type="model output",q05n=qq[1,],q95n=qq[3,],
                                                q05=pr$lower95,q95=pr$upper95)
                            } else
                            {
                              df  <- data.frame(y=qq[2,],type="model output",q05n=qq[1,],q95n=qq[3,],
                                                q05=pr$lower95,q95=pr$upper95)
                            }
                          } else if (CI == "GP")
                          {
                            nugget <- mvrnorm(n=100,pr$mean,diag(var,100))
                            qq <- apply(nugget,2,quantile,c(0.05,0.5,0.95))
                            df  <- data.frame(y=qq[2,],type="model output",q05=pr$lower95,
                                              q95=pr$upper95, fill="CI 90% GP")
                          } else
                          {
                            warning("The argument for the credibility interval is not valid and no credibility interval
                                                        will be displayed")
                            df <- data.frame(y=pr$mean+rnorm(nrow(X),0,sqrt(var)),
                                             type="model output")
                          }
                          self$cov.gp <- pr$cov
                          self$yc <- pr$mean
                          return(df)
                        })
                        )


model2.class$set("public","surrogate",
                 function()
                 {
                   if (self$case == "1")
                   {
                     ## Dim is the dimension of H*Theta
                     Dim       <- self$opt.emul$p+self$d
                     ## Creation of the maximin LHS
                     doe       <- lhsDesign(self$opt.emul$n.emul,Dim)$design
                     doe       <- maximinSA_LHS(doe)$design
                     ## Get the boundaries of X
                     binf.X    <- apply(self$X,2,min)
                     bsup.X    <- apply(self$X,2,max)
                     ## Going back to the original space H and Theta
                     doe.X     <- unscale(doe[,c(1:self$d)],binf.X,bsup.X)
                     doe.theta <- unscale(doe[,c((self$d+1):Dim)],self$opt.emul$binf,
                                            self$opt.emul$bsup)
                     ## Generate the doe
                     self$doe  <- cbind(doe.X,doe.theta)
                     ## Compute the output of the code for the doe
                     self$z <- NULL
                     for (i in 1:self$opt.emul$n.emul)
                     {
                       covariates <- as.matrix(self$doe[i,1:(Dim-self$opt.emul$p)])
                       dim(covariates) <- c(1,self$d)
                       self$z <- c(self$z,self$code(covariates,self$doe[i,(Dim-self$opt.emul$p+1):Dim]))
                     }
                     ## Create the Gaussian Process corresponding
                     GP <- km(formula =~1, design=self$doe, response = self$z,covtype = self$opt.gp$type)
                     return(GP)

                   } else if (self$case == "2")
                   {
                     self$doe <- self$opt.gp$DOE
                     Dim <- ncol(self$doe)
                     self$p <- Dim - self$d
                     ### Generating the response
                     self$z <- NULL
                     for (i in 1:nrow(self$doe))
                     {
                       covariates <- as.matrix(self$doe[i,1:(self$d)])
                       dim(covariates) <- c(1,self$d)
                       self$z <- c(self$z,self$code(covariates,self$doe[i,(Dim-self$p+1):Dim]))
                     }
                     ## Create the Gaussian Process corresponding
                     GP <- km(formula =~1, design=self$doe, response = self$z,covtype = self$opt.gp$type)
                     return(GP)
                   } else if (self$case == "3")
                   {
                     self$p <- ncol(self$opt.sim$DOEsim)-self$d
                     ## GP from design of expiriments and code outputs
                     GP <- km(formula =~1, design=as.data.frame(self$opt.sim$DOEsim),
                              response = self$opt.sim$Ysim,covtype = self$opt.gp$type)
                     return(GP)
                   }
                 })


## Likelihood
model2.class$set("public","likelihood",
                 function(theta,var)
                 {
                   self$m.exp <- self$yc
                   self$V.exp <- var*diag(self$n) + self$cov.gp
                   return(-self$n/2*log(2*pi)-1/2*log(det(self$V.exp))
                          -0.5*t(self$Yexp-self$m.exp)%*%solve(self$V.exp)%*%(self$Yexp-self$m.exp))
                 })


## plot function
model2.class$set("public","plot",
                 function(x,CI="all",...)
                 {
                   ### Plot generates a ggplot object
                   if (is.matrix(x)){stop("please enter a correct x to plot your model")}
                   if (length(x)!= self$n){stop(paste("please enter a correct vector x of size",
                                                      self$n,sep=" "))}
                   df <- data.frame(y=self$Yexp,type="exp")
                   if (!is.null(self$theta) & !is.null(self$var)){
                     df2 <- cbind(self$model.fun(self$theta,self$var,self$X,CI),x=x)
                     if (!is.null(CI))
                     {
                       if(CI == "err" | CI == "all")
                       {
                         df.noise <- cbind(df,q05=df2$q05n,q95=df2$q95n,fill="CI 90% noise",x=x)
                         df2.noise <- data.frame(y=df2$y,type=df2$type,q05=df2$q05n,
                                                 q95=df2$q95n,fill="CI 90% noise",x=x)
                       }
                       if (CI == "GP" | CI =="all")
                       {
                         df.gp <- cbind(df,q05=df2$q05,q95=df2$q95,fill="CI 90% GP",x=x)
                         df2.gp <- data.frame(y=df2$y,type=df2$type,q05=df2$q05,q95=df2$q95,fill="CI 90% GP",x=x)
                       }
                     }
                     if (is.null(CI))
                     {
                       df  <- cbind(df,x=x)
                       df  <- rbind(df,df2)
                       p <- ggplot(df) + geom_line(mapping = aes(x=x,y=y,color=type))+
                         xlab("")+ylab("")+theme_light()+self$gglegend()+
                         scale_color_manual(values=c("red", "#000000"))
                     } else if (CI == "err")
                     {
                       df   <- rbind(df.noise,df2.noise)
                       p <- ggplot(df)+geom_line(mapping = aes(x=x,y=y, color=type))+
                         xlab("")+ylab("")+theme_light()+self$gglegend()+
                         geom_ribbon(mapping = aes(x=x,ymin=q05,ymax=q95,fill="CI 90% noise"),
                                     alpha=0.4)+scale_fill_manual(values = adjustcolor("skyblue3"))+
                         scale_color_manual(values=c("red", "#000000"))
                     } else if (CI == "GP")
                     {
                       df   <- rbind(df.gp,df2.gp)
                       p <- ggplot(df)+geom_line(mapping = aes(x=x,y=y, color=type))+
                         xlab("")+ylab("")+theme_light()+self$gglegend()+
                         geom_ribbon(mapping = aes(x=x,ymin=q05,ymax=q95,fill="CI 90% GP"),
                                     alpha=0.4)+scale_fill_manual(values = adjustcolor("grey12"))+
                         scale_color_manual(values=c("red", "#000000"))
                     } else if (CI == "all")
                     {
                       df  <- rbind(df.gp,df2.gp)
                       df2 <- rbind(df.noise,df2.noise)
                       p <- ggplot(df)+
                         xlab("")+ylab("")+theme_light()+self$gglegend()+
                         geom_ribbon(mapping = aes(x=x,ymin=q05,ymax=q95,fill="CI 90% noise"),
                                     alpha=0.3)+
                         geom_ribbon(data=df2,mapping = aes(x=x,ymin=q05,ymax=q95,fill="CI 90% GP"),
                                     alpha=0.8,show.legend = FALSE)+
                         geom_line(mapping = aes(x=x,y=y, color=type))+
                         scale_fill_manual(name = NULL,values = adjustcolor(c("skyblue3", "grey12"),
                         alpha.f = 0.3))+
                         guides(fill = guide_legend(override.aes = list(alpha = c(0.8,0.3))))+
                         scale_color_manual(values=c("red", "#000000"))
                     } else
                     {
                       warning("The argument for the credibility interval is not valid and no credibility interval
                                                        will be displayed")
                       df  <- cbind(df,x=x)
                       df  <- rbind(df,df2)
                       p <- ggplot(df) + geom_line(mapping = aes(x=x,y=y,color=type))+
                         xlab("")+ylab("")+theme_light()+self$gglegend()+
                         scale_color_manual(values=c("red", "#000000"))
                     }
                   } else
                   {
                     p <- ggplot(df) + geom_line(mapping = aes(x=x,y=y,color=type))+
                       xlab("")+ylab("")+theme_light()+self$gglegend()+
                       scale_color_manual(values=c("red"))
                   }
                   return(p)
                 })


### Print function
model2.class$set("public","print",
                 function()
                 {
                   cat("Call:\n")
                   print(self$model)
                   cat("\n")
                   cat("With the function:\n")
                   print(self$code)
                   cat("\n")
                   cat("A surrogate had been set up:")
                   print(self$GP)
                   cat("\n")
                   cat("No discrepancy is added\n")
                 }
)


################################## Model 4 definition ####################################

## model 4 main functions
model4.class <- R6Class(classname = "model4.class",
                        inherit = model2.class,
                        public=list(
                          funC = NULL,
                          predTemp = NULL,
                          disc = NULL,
                          initialize=function(code=NA, X=NA, Yexp=NA, model=NA,opt.gp=NA,opt.emul=NA,
                                              opt.disc=list(kernel.type=NULL),opt.sim=NA)
                          {
                            super$initialize(code, X, Yexp, model,opt.gp, opt.emul,opt.sim)
                            self$opt.disc  <- list(kernel.type=opt.disc$kernel.type)
                            if (is.null(self$opt.disc$kernel.type)==TRUE)
                            {
                              self$opt.disc$kernel.type="gauss"
                            }
                            self$funC     <- super$fun
                            self$predTemp <- super$pred
                          },
                          discrepancy = function(theta,thetaD,var,X=self$X)
                          {
                            y   <- self$funC(theta,var)$y
                            z   <- self$Yexp - y
                            Cov <- kernel.fun(X,thetaD[1],thetaD[2],self$opt.disc$kernel.type)
                            if (is.null(dim(X)) && length(X)==1)
                            {} else
                            {
                              p <- eigen(Cov)$vectors
                              e <- eigen(Cov)$values
                              if (all(e>0)){} else
                              {
                                e[which(e<0)] <- 1e-4
                              }
                              d <- diag(e)
                              if (nrow(p) == 1 & ncol(p) == 1)
                              {
                                Cov <- as.numeric(p)^2*d
                              } else
                              {
                                Cov <- t(p)%*%d%*%p
                              }
                            }
                            if (is.null(dim(X))){long <- length(X)}else{
                              long <- dim(X)[1]}
                            if (long==1)
                            {
                              if (nrow(p) == 1 & ncol(p) == 1)
                              {
                                bias <- rnorm(1,0,sqrt(Cov))
                              } else
                              {
                                bias <- rnorm(n=self$n,0,sqrt(Cov))
                              }
                            } else
                            {
                              bias <- mvrnorm(n=self$n,rep(0,long),Cov)
                              bias <- apply(bias,1,mean)
                            }
                            return(list(bias=bias,cov=Cov))
                          },
                          pred = function(theta,thetaD,var,x.new)
                          {
                            self$disc <- self$discrepancy(theta,thetaD,var,x.new)
                            foo <- self$predTemp(theta,var,x.new)
                            y <- foo$y
                            yc  <- foo$yc
                            return(list(y=self$disc$bias+y,cov=self$disc$cov,yc=yc))
                          },
                          fun = function(theta,thetaD,var,X=self$X)
                          {
                            foo <- self$funC(theta,var)
                            self$disc <- self$discrepancy(theta,thetaD,var,X)
                            y <- foo$y
                            Cov.GP <- foo$Cov.GP
                            yc <- foo$yc
                            lower <- foo$lower
                            upper <- foo$upper
                            return(list(y=self$disc$bias+y,Cov.D=self$disc$cov,yc=yc,
                                        lower=lower,upper=upper,Cov.GP=Cov.GP))
                          })
)


model4.class$set("public","likelihood",
                 function(theta,thetaD,var)
                 {
                   temp <- self$fun(as.vector(theta),thetaD,var)
                   self$m.exp <- temp$yc
                   self$V.exp <- var*diag(self$n) + temp$Cov.GP +temp$Cov.D
                   # return(1/((2*pi)^(self$n/2)*det(self$V.exp)^(1/2))*exp(-1/2*t(self$Yexp-self$m.exp)%*%
                   #                                                invMat(self$V.exp)%*%(self$Yexp-self$m.exp)))
                   return(-0.5*t(self$Yexp-self$m.exp)%*%solve(self$V.exp)%*%(self$Yexp-self$m.exp))
                 })



model4.class$set("public","plot",
                 function(theta,thetaD,var,select.X=NULL,CI=c("GP","err"),points=FALSE)
                 {
                   if (length(theta)!=self$p)
                   {stop('You have given the wrong number of parameter')}
                   if(self$d>1 & is.null(select.X))
                   {stop('Graphic representation is not available in dimension >1')}
                   if(is.null(select.X))
                   {
                     res <- self$fun(theta,thetaD,var,self$X)
                   } else{
                     res <- self$fun(theta,thetaD,var,select.X)
                     }
                   Xplot <- self$X
                   if(is.null(select.X)==FALSE){Xplot <- select.X}
                   binf <- min(Xplot)
                   bsup <- max(Xplot)
                   funCpp <- function(theta,thetaD,var)
                   {
                     return(self$fun(theta,thetaD,var)$y)
                   }
                   if (length(CI)==1)
                   {
                     if (CI=="GP")
                     {
                       gg.data <- data.frame(y=res$yc,x=seq(binf,bsup,length.out=length(res$yc)),
                                             lower=res$lower,upper=res$upper,type="Gaussian process",
                                             fill="CI 90% GP")
                       gg.data.exp <- data.frame(y=self$Yexp,x=seq(binf,bsup,length.out=length(res$yc)),
                                                 lower=res$lower,
                                                 upper=res$upper,type="Experiment",
                                                 fill="CI 90% GP")
                       gg.data <- rbind(gg.data,gg.data.exp)
                       if (is.null(self$code))
                       {
                         if (is.null(select.X))
                         {
                           gg.points <- data.frame(x=self$DOEsim[,1],y=self$Ysim)
                         } else
                         {
                           gg.points <- data.frame(x=self$Xplot,y=self$Ysim)
                         }
                       }else
                       {
                         if (is.null(select.X))
                         {
                           gg.points <- data.frame(x=self$D[,1],y=self$Yc)
                         } else
                         {
                           gg.points <- data.frame(x=self$Xplot,y=self$Yc)
                         }
                       }
                       p <- ggplot(gg.data)+ geom_ribbon(aes(ymin=lower,ymax=upper,x=x,fill=fill),alpha=0.3)+
                         geom_line(aes(y=y,x=x,col=type))+
                         theme_light()+
                         ylab("")+xlab("")+
                         scale_fill_manual("",values=c("grey12"))+
                         theme(legend.position=c(0.65,0.86),
                               legend.text=element_text(size = '12'),
                               legend.title=element_blank(),
                               legend.key=element_rect(colour=NA),
                               axis.text=element_text(size=20))
                       if (points==FALSE)
                       {
                         return(p)
                       } else
                       {
                         return(p+geom_jitter(data=gg.points,aes(x=x,y=y)))
                       }
                     } else {
                       if (CI=="err")
                         yres <- resCppD(funCpp,theta,thetaD,var)
                         qqerr <- apply(yres,1,quantile,c(0.05,0.95))
                         gg.data <- data.frame(y=res$y,x=seq(binf,bsup,length.out=length(res$yc)),
                                               type="Model output")
                         gg.data.n <- data.frame(x=seq(binf,bsup,length.out=length(res$yc)),
                                                 ymin=qqerr[1,],ymax=qqerr[2,],
                                               type="CI 90% discrepancy + noise")
                         gg.data.exp  <- data.frame(y=self$Yexp,x=seq(binf,bsup,length.out=length(res$yc)),
                                                  type="Experiments")
                         gg.data <- rbind(gg.data,gg.data.exp)
                         p <- ggplot(gg.data) +
                            geom_ribbon(data = gg.data.n,aes(x=x,ymin=ymin,ymax=ymax,fill=type),alpha=0.8)+
                            geom_line(aes(y=y,x=x,col=type)) +
                            theme_light() + scale_fill_manual(values = "skyblue3")+
                            ylab("")+xlab("") +
                            theme(legend.position=c(0.65,0.86),
                               legend.text=element_text(size = '12'),
                               legend.title=element_blank(),
                               legend.key=element_rect(colour=NA),
                               axis.text=element_text(size=20))
                       return(p)
                     }
                   } else {
                     if (length(CI)==2)
                     {
                       if ((CI[1]=="GP" & CI[2]=="err") | (CI[2]=="GP" & CI[1]=="err"))
                       {
                         yres <- resCppD(funCpp,theta,thetaD,var)
                         qqerr <- apply(yres,1,quantile,c(0.05,0.95))
                         gg.data <- data.frame(y=res$yc,x=seq(binf,bsup,length.out=length(res$yc)),
                                               lower=res$lower,upper=res$upper,type="Gaussian process",
                                               fill="CI 90% GP")
                         gg.data.exp <- data.frame(y=self$Yexp,x=seq(binf,bsup,length.out=length(res$yc)),
                                                   lower=res$lower,
                                                   upper=res$upper,type="Experiment",
                                                   fill="CI 90% GP")
                         gg.data.n <- data.frame(x=seq(binf,bsup,length.out=length(res$yc)),
                                                 ymin=qqerr[1,],ymax=qqerr[2,],
                                                 type="CI 90% discrepancy + noise")
                         gg.data <- rbind(gg.data,gg.data.exp)
                         if (is.null(self$code))
                         {
                           if (is.null(select.X))
                           {
                             gg.points <- data.frame(x=self$DOEsim[,1],y=self$Ysim)
                           } else
                           {
                             gg.points <- data.frame(x=self$Xplot,y=self$Ysim)
                           }
                         }else
                         {
                           if (is.null(select.X))
                           {
                             gg.points <- data.frame(x=self$D[,1],y=self$Yc)
                           } else
                           {
                             gg.points <- data.frame(x=self$Xplot,y=self$Yc)
                           }
                         }
                         p <- ggplot(gg.data)+ geom_ribbon(aes(ymin=lower,ymax=upper,x=x,fill=fill),alpha=0.3)+
                           geom_ribbon(data = gg.data.n,aes(x=x,ymin=ymin,ymax=ymax,fill=type),
                                       alpha=0.8,show.legend = FALSE)+
                           geom_line(aes(y=y,x=x,col=type))+
                           theme_light()+
                           ylab("")+xlab("")+
                           scale_fill_manual(name = NULL,values = adjustcolor(c("skyblue3", "grey12"),
                                                                              alpha.f = 0.3))+
                           guides(fill = guide_legend(override.aes = list(alpha = c(0.8,0.3)))) +
                           theme(legend.position=c(0.65,0.86),
                                 legend.text=element_text(size = '12'),
                                 legend.title=element_blank(),
                                 legend.key=element_rect(colour=NA),
                                 axis.text=element_text(size=20))
                         if (points==FALSE)
                         {
                           return(p)
                         } else
                         {
                           return(p+geom_jitter(data=gg.points,aes(x=x,y=y)))
                         }
                       } else
                       {print("Enter the right value for the CI option")}
                     } else
                     {print("The length of the CI option is too long")}
                   }
                 })

model4.class$set("public","print",
                 function()
                 {
                   cat("Call:\n")
                   print(self$model)
                   cat("\n")
                   cat("With the function:\n")
                   print(self$code)
                   cat("\n")
                   cat("A surrogate had been set up:")
                   print(self$GP)
                   cat("\n")
                   cat("A discrepancy is added:\n")
                   # cat(paste("Mean of the bias:",round(mean(self$disc$bias),5),"\n",sep=" "))
                   # cat(paste("Covariance of the bias:",round(mean(self$disc$cov),5),"\n",sep=" "))
                   cat("Chosen kernel:", self$opt.disc$kernel.type)
                 }
)

