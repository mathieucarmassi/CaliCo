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
#' @field opt.pg a list of parameter for the surrogate (default NULL) \itemize{
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
                   initialize = function(code=NA,X=NA,Yexp=NA,model=NA)
                   {
                     self$code  <- code
                     self$X     <- X
                     self$Yexp  <- Yexp
                     self$n     <- length(Yexp)
                     self$model <- model
                     if (is.matrix(X)) {self$d <- ncol(X)} else{self$d <-1}
                     private$checkModels()
                     private$checkCode()
                     private$checkOptions()
                   },
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
                    test(c("p","n.emul","binf","bsup"),self$opt.emul)
                    test(c("type","DOE"),self$opt.pg)
                    test(c("Ysim","DOEsim"),opt.sim)
                  } else
                  {
                    if (self$model %in% c("model3","model4")){test(c("Kernel.type"),opt.disc)}
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


model1.class <- R6Class(classname = "model1.class",
                        inherit = model.class,
                        public=list(
                        m.exp = NULL, ## Mean for the likelihood
                        V.exp = NULL, ## Variance for the likelihood
                        l     = NULL, ## Length of the new data for prediction
                        theta = NULL, ## parameter vector
                        var   = NULL, ## variance of the measurement error
                        initialize=function(code=NA, X=NA, Yexp=NA, model=NA)
                        {
                          ### Initialize from model.class
                          super$initialize(code, X, Yexp, model)
                        },
                        model.fun = function(theta,var,CI=FALSE)
                        {
                          ### Function that generates the output of the model. If CI=TRUE, it computes the credibility
                          ### intervals of the white Gaussian noise
                          if (CI==TRUE)
                          {
                            y <- matrix(nr=100,nc=self$n)
                            for(i in 1:100){y[i,] <- self$code(self$X,theta)+rnorm(self$n,0,sqrt(var))}
                            qq <- apply(y,2,quantile,c(0.05,0.5,0.95))
                            df <- data.frame(y=qq[2,],type="code output",q05=qq[1,],q95=qq[3,],fill="CI 90% noise")
                          }else
                          {
                            y  <- self$code(self$X,theta)+rnorm(self$n,0,sqrt(var))
                            df <- data.frame(y=y,type="code output")
                          }
                          return(df)
                        },
                        prediction.fun = function(theta,var,x.new)
                        {
                          ### Prediction function is the function to use when applying on a new data set
                          if (is.matrix(x.new)){l <- nrow(x.new)} else{l <- length(x.new)}
                          y  <- self$code(x.new,theta)+rnorm(l,0,sqrt(var))
                          df <- data.frame(pr=y,type="predicted")
                          return(df)
                        },
                        likelihood = function(theta,var)
                        {
                          ### Likelihood
                          self$m.exp <- self$code(self$X,as.vector(theta))
                          self$V.exp <- var*diag(self$n)
                          return(-0.5*t(self$Yexp-self$m.exp)%*%solve(self$V.exp)%*%(self$Yexp-self$m.exp))
                        },
                        print = function(...)
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
                        },
                        plot = function(x,theta=NULL,var=NULL,CI="noise",...)
                        {
                          ### Plot generates a ggplot object
                          if (is.matrix(x)){stop("please enter a correct x to plot your model")}
                          if (length(x)!= self$n){stop(paste("please enter a correct vector x of size",
                                                             self$n,sep=" "))}
                          df <- data.frame(y=self$Yexp,type="exp")
                          if (!is.null(theta) & !is.null(var)) {
                            if (is.null(CI)){CI=""}
                            if (CI=="noise")
                            {
                              df2 <- self$model.fun(theta,var,CI=TRUE)
                              df2 <- cbind(df2,x=x)
                              df  <- cbind(df,q05=df2$q05,q95=df2$q95,fill="CI 90% noise",x=x)
                              df  <- rbind(df,df2)
                              p <- ggplot(df) + geom_line(mapping = aes(x=x,y=y,color=type))+
                                xlab("")+ylab("")+theme_light()+self$gglegend()+
                                geom_ribbon(mapping = aes(x=x,ymin=q05,ymax=q95,fill=fill),alpha=0.4)+
                                scale_fill_manual(values = "grey70")
                            } else {
                              df2 <- self$model.fun(theta,var,CI=FALSE)
                              df2 <- cbind(df2,x=x)
                              df  <- cbind(df,x=x)
                              df  <- rbind(df,df2)
                              p <- ggplot(df) + geom_line(mapping = aes(x=x,y=y,color=type))+
                                xlab("")+ylab("")+theme_light()+self$gglegend()
                            }
                          }
                          return(p)
                        }
                        )
                        )



model3.class <- R6Class(classname = "model3.class",
                        inherit = model1.class,
                        public=list(
                          model1.fun        = NULL, ## model function from model1
                          model1.prediction = NULL, ## prediction function from model1
                          opt.disc          = NULL, ## discrepancy options
                          disc              = NULL, ## discrepancy field
                          initialize=function(code=NA, X=NA, Yexp=NA, model=NA,opt.disc=list(type.kernel=NULL))
                          {
                            ## Check if the opt.emul option is filled if it is not a gaussian kernel is picked
                            if (is.null(opt.disc$type.kernel))
                            {
                              warning("default value is selected. The discrepancy will have a Gaussian covariance
                                      structure")
                              self$opt.disc$kernel.type="gauss"
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
                          ## Define method that generates a discrepancy
                          discrepancy = function(theta,thetaD,var,X=self$X)
                          {
                            y   <- self$model(theta,var)$yfun
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
                                biais <- rnorm(1,0,sqrt(Cov))
                              } else
                              {
                                biais <- rnorm(n=self$n,0,sqrt(Cov))
                              }
                            } else
                            {
                              biais <- mvrnorm(100,rep(0,long),Cov)
                              dim(biais) <- c(long,100)
                              biais <- apply(biais,1,mean)
                            }
                            return(list(biais=biais,cov=Cov))
                          },
                          model.fun = function(theta,thetaD,var,X=self$X)
                          {
                            self$disc <- self$discrepancy(theta,thetaD,var,X)
                            foo <- self$model1.fun(theta,var)
                            y <- foo$y
                            yc  <- foo$yc
                            return(list(y=self$disc$biais+y,cov=self$disc$cov,yc=yc))
                          },
                          prediction.fun = function(theta,thetaD,var,x.new)
                          {
                            self$disc <- self$discrepancy(theta,thetaD,var,x.new)
                            foo <- self$model1.fun(theta,var,x.new)
                            y <- foo$y
                            yc  <- foo$yc
                            return(list(y=self$disc$biais+y,cov=self$disc$cov,yc=yc))
                          }
                          )
)


model3.class$set("public","plot",
                 function(theta,thetaD,var,select.X=NULL,CI=TRUE)
                 {
                   if(is.null(dim(self$X))==FALSE & is.null(select.X))
                   {stop('Graphic representation is not available in dimension >1')}
                   res <- self$fun(theta,thetaD,var)
                   Xplot <- self$X
                   if(is.null(select.X)==FALSE){Xplot <- select.X}
                   binf <- min(Xplot)
                   bsup <- max(Xplot)
                   if (CI==FALSE)
                   {
                     gg.data <- data.frame(y=res$yc,x=seq(binf,bsup,length.out=length(res$yc)),type="Model Output")
                     gg.data.exp  <- data.frame(y=self$Yexp,x=seq(binf,bsup,length.out=length(res$yc)),
                                                type="Experiments")
                     gg.data <- rbind(gg.data,gg.data.exp)
                     p <- ggplot(gg.data) +
                       geom_line(aes(y=y,x=x,col=type)) +
                       theme_light() +
                       ylab("")+xlab("") +
                       theme(legend.position=c(0.65,0.86),
                             legend.text=element_text(size = '12'),
                             legend.title=element_blank(),
                             legend.key=element_rect(colour=NA),
                             axis.text=element_text(size=20))
                   } else
                   {
                     funCpp <- function(theta,thetaD,var)
                     {
                       return(self$fun(theta,thetaD,var)$y)
                     }
                     yres <- resCppD(funCpp,theta,thetaD,var)
                     qqres <- apply(yres,1,quantile,c(0.05,0.95))
                     gg.data <- data.frame(y=res$yc,x=seq(binf,bsup,length.out=length(res$yc)),type="Model Output")
                     gg.data.n <- data.frame(x=seq(binf,bsup,length.out=length(res$yc)),ymin=qqres[1,],ymax=qqres[2,],
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
                   }
                   return(p)
                 })

model3.class$set("public","print",
                 function()
                 {
                   cat("Call:\n")
                   print(self$model)
                   cat("\n")
                   cat("With the function:\n")
                   print(self$code)
                   cat("\n")
                   cat("No surrogate is selected")
                   cat("\n\n")
                   cat("A discrepancy is added:\n")
                   cat(paste("Mean of the biais:",round(mean(self$disc$biais),5),"\n",sep=" "))
                   cat(paste("Covariance of the biais:",round(mean(self$disc$cov),3),"\n",sep=" "))
                   cat(paste("Kernel chossen: ",self$opt.disc$kernel.type,sep=""))
                 }
)



model3.class$set("public","likelihood",
                 function(theta,thetaD,var)
                 {
                   self$m.exp <- self$code(self$X,as.vector(theta))
                   temp <- self$fun(as.vector(theta),thetaD,var)
                   self$V.exp <- var*diag(self$n) + temp$cov
                   # return(1/((2*pi)^(self$n/2)*det(self$V.exp)^(1/2))*exp(-1/2*t(self$Yexp-self$m.exp)%*%
                   #                                              invMat(self$V.exp)%*%(self$Yexp-self$m.exp)))
                   return(-0.5*t(self$Yexp-self$m.exp)%*%solve(self$V.exp)%*%(self$Yexp-self$m.exp))
                 })



model2.class <- R6Class(classname = "model2.class",
                        inherit = model.class,
                        public = list(
                          n.emul   = NULL,
                          opt.emul = NULL,
                          opt.sim  = NULL,
                          Ysim     = NULL,
                          DOEsim   = NULL,
                          GP       = NULL,
                          X        = NULL,
                          DOE      = NULL,
                          design   = NULL,
                          p        = NULL,
                          type     = NULL,
                          D        = NULL,
                          m.exp    = NULL,
                          V.exp    = NULL,
                        initialize = function(code=NA, X=NA, Yexp=NA, model=NA,opt.pg=list(type=NA,DOE=NA),
                                              opt.emul=list(p=NA,n.emul=NA,binf=NA,bsup=NA),
                                              opt.sim=list(Ysim=NA,DOEsim=NA))
                        {
                          super$initialize(code, X, Yexp, model,opt.emul)
                          self$opt.emul <- opt.emul
                          self$X        <- X
                          self$binf     <- opt.emul$binf
                          self$bsup     <- opt.emul$bsup
                          self$n.emul   <- opt.emul$n.emul
                          self$p        <- opt.emul$p
                          self$type     <- opt.pg$type
                          self$DOE      <- opt.pg$DOE
                          self$opt.sim  <- opt.sim
                          self$Ysim     <- opt.sim$Ysim
                          self$DOEsim   <- opt.sim$DOEsim
                          private$checkOptions()
                          if (is.null(self$code))
                          {
                            if (is.null(opt.sim))
                            {
                              print("The numerical code is desabled, please fill the opt.sim option")
                            }
                          }
                          foo           <- self$surrogate()
                          self$GP       <- foo$GP
                          self$design   <- foo$doe
                          print("The surrogate has been set up, you can now use the function")
                        },
                        surrogate = function()
                        {
                          if (is.null(self$code))
                          {
                            GP <- km(formula =~1, design=as.data.frame(self$DOEsim),
                                     response = self$Ysim,covtype = self$type)
                            return(list(GP=GP))
                          } else
                          {
                            Xcr <- scale(self$X)
                            V   <- attr(Xcr,"scaled:scale")
                            M   <- attr(Xcr,"scaled:center")
                            Dim <- self$p+self$d
                            if (is.null(self$DOE)==FALSE)
                            {
                              self$n.emul <- dim(self$DOE)[1]
                              self$D <- self$DOE
                              ### Generating the response
                              z <- matrix(nr=self$n.emul,nc=1)
                              for (i in 1:self$n.emul)
                              {
                                covariates <- as.matrix(self$D[i,1:(Dim-self$p)])
                                dim(covariates) <- c(1,self$d)
                                z[i] <- self$code(covariates,self$D[i,(Dim-self$p+1):Dim])
                              }
                              #z <- self$code(D[,1:(Dim-self$p)],D[,(Dim-self$p+1):Dim])
                              ### Converting D as a data.frame for the km function
                              self$D <- as.data.frame(self$D)
                              ### Creation of the Gaussian Process with estimation of hyperpameters
                              GP <- km(formula =~1, design=self$D, response = z,covtype = self$type)
                              return(list(GP=GP,doe=self$D))
                            }else
                            {
                              doe <- lhsDesign(self$n.emul,Dim)$design
                              doe <- maximinSA_LHS(doe)
                              doe <- doe$design
                              ### Getting back the value of the parameter generated by the DOE
                              binf.X <- apply(Xcr,2,min)
                              bsup.X <- apply(Xcr,2,max)
                              DOEtemp <- unscale(doe[,1:self$d],binf.X,bsup.X)
                              if (self$d==1)
                              {
                                for (i in 1:self$n.emul)
                                {
                                  DOEtemp[i] <- DOEtemp[i]*V+M
                                }
                              } else {
                                for (i in 1:self$n.emul)
                                {
                                  DOEtemp[i,] <- DOEtemp[i,]*V+M
                                }
                              }
                              if (length(self$binf)!=self$p | length(self$bsup)!=self$p)
                              {stop('Mismatch between the size of the upper, lower bounds and the parameter vector')}
                              doeParam <- unscale(doe[,(Dim-self$p+1):Dim],self$binf,self$bsup)
                              ### Matrix D contains the final value for the DOE
                              self$D <- cbind(DOEtemp,doeParam)
                              self$DOE <- self$D
                              ### Generating the response
                              z <- matrix(nr=self$n.emul,nc=1)
                              for (i in 1:self$n.emul)
                              {
                                covariates <- as.matrix(self$D[i,1:(Dim-self$p)])
                                dim(covariates) <- c(1,self$d)
                                z[i] <- self$code(covariates,self$D[i,(Dim-self$p+1):Dim])
                              }
                              #z <- self$code(D[,1:(Dim-self$p)],D[,(Dim-self$p+1):Dim])
                              ### Converting D as a data.frame for the km function
                              self$D <- as.data.frame(self$D)
                              #colnames(D) <- c("V1","V2","V3","V4","V5","V6")
                              ### Creation of the Gaussian Process with estimation of hyperpameters
                              GP <- km(formula =~1, design=self$D, response = z,covtype = self$type)
                              return(list(GP=GP,doe=self$D))
                            }
                          }
                        },
                        fun = function(theta,var)
                        {
                          if(is.null(dim(self$X)) && length(self$X)==1 && self$X==0)
                          {
                            if(self$p==1)
                            {
                              Xnew <- rep(theta,self$n)
                            } else
                            {
                              Xnew <- matrix(rep(theta,rep(self$n,self$p)),nr=self$n,nc=self$p)
                            }
                          } else
                          {
                            if(self$p==1)
                            {
                              Xnew <- cbind(self$X,rep(theta,self$n))
                            } else
                            {
                              Xtemp <- matrix(rep(theta,rep(self$n,self$p)),nr=self$n,nc=self$p)
                              Xnew  <- cbind(self$X,Xtemp)
                            }
                          }
                          Xnew <- as.data.frame(Xnew)
                          names(Xnew) <- c("DOE","doeParam")
                          pr <- predict(self$GP,newdata=as.data.frame(Xnew),type="UK",
                                        cov.compute=TRUE,interval="confidence",checkNames=FALSE)
                          err <- rnorm(n=self$n,mean = 0,sd=sqrt(var))
                          return(list(y=pr$mean+err,Cov.GP=pr$cov,yc=pr$mean,lower=pr$lower95,upper=pr$upper95))
                        },
                        pred = function(theta,var,x.new)
                        {
                          if (is.matrix(x.new)){l <- nrow(x.new)} else{l <- length(x.new)}
                          if(self$p==1)
                          {
                            Xnew <- cbind(x.new,rep(theta,l))
                          } else
                          {
                            Xtemp <- matrix(rep(theta,rep(l,self$p)),nr=l,nc=self$p)
                            Xnew  <- cbind(x.new,Xtemp)
                          }
                          Xnew <- as.data.frame(Xnew)
                          #names(Xnew) <- c("DOE","doeParam")
                          pr <- predict(self$GP,newdata=as.data.frame(Xnew),type="UK",
                                        cov.compute=TRUE,interval="confidence",checkNames=FALSE)
                          if (length(pr$mean) == 1)
                          {
                            err <- rnorm(n=1,mean = 0,sd=sqrt(var))
                          } else
                          {
                            err <- rnorm(n=self$n,mean = 0,sd=sqrt(var))
                          }
                          return(list(y=pr$mean+err,Cov.GP=pr$cov,yc=pr$mean,lower=pr$lower95,upper=pr$upper95))
                        })
                        )


model2.class$set("public","likelihood",
                 function(theta,var)
                 {
                   temp <- self$fun(as.vector(theta),var)
                   self$m.exp <- temp$yc
                   if (length(theta)!=self$p)
                   {stop('You have given the wrong number of parameter')}
                   self$V.exp <- var*diag(self$n) + temp$Cov.GP
                   # return(1/((2*pi)^(self$n/2)*det(self$V.exp)^(1/2))*exp(-1/2*t(self$Yexp-self$m.exp)%*%
                   #                                                invMat(self$V.exp)%*%(self$Yexp-self$m.exp)))
                   return(-0.5*t(self$Yexp-self$m.exp)%*%solve(self$V.exp)%*%(self$Yexp-self$m.exp))
                 })


model2.class$set("public","plot",
                 function(theta, var, select.X=NULL, CI=c("GP","err"), points=FALSE)
                 {
                   if (length(theta)!=self$p)
                   {stop('You have given the wrong number of parameter')}
                   if(self$d>1 & is.null(select.X))
                   {stop('Graphic representation is not available in dimension >1')}
                   res <- self$fun(theta,var)
                   funCpp <- function(theta,var)
                   {
                     return(self$fun(theta,var)$y)
                   }
                   Xplot <- self$X
                   if(is.null(select.X)==FALSE){
                     Xplot <- select.X
                   }
                   binf <- min(Xplot)
                   bsup <- max(Xplot)
                   if (length(CI)==1)
                   {
                     if (CI=="GP")
                     {
                       gg.data <- data.frame(y=res$yc,x=seq(binf,bsup,length.out=length(res$yc)),
                                             lower=res$lower,upper=res$upper,type="Gaussian process",
                                             fill="CI 90% GP")
                       gg.data.exp <- data.frame(y=self$Yexp,x=seq(binf,bsup,length.out=length(res$yc)),
                                                 lower=res$lower,upper=res$upper,type="Experiment",
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
                          yres <- resCpp(funCpp,theta,var)
                          qqerr <- apply(yres,1,quantile,c(0.05,0.95))
                          gg.data <- data.frame(y=res$y,x=seq(binf,bsup,length.out=length(res$yc)),
                                                type="Model output")
                          gg.data.n <- data.frame(x=seq(binf,bsup,length.out=length(res$yc)),
                                                  ymin=qqerr[1,],ymax=qqerr[2,],type="CI 90% noise")
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
                         yres <- resCpp(funCpp,theta,var)
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
                                                 type="CI 90% noise")
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
                           scale_fill_manual(name = NULL,values = adjustcolor(c("grey12", "skyblue3"),
                                                                              alpha.f = 0.3))+
                           guides(fill = guide_legend(override.aes = list(alpha = c(0.3,0.8)))) +
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



model4.class <- R6Class(classname = "model4.class",
                        inherit = model2.class,
                        public=list(
                          funC = NULL,
                          predTemp = NULL,
                          disc = NULL,
                          initialize=function(code=NA, X=NA, Yexp=NA, model=NA,opt.pg=NA,opt.emul=NA,
                                              opt.disc=list(kernel.type=NULL),opt.sim=NA)
                          {
                            super$initialize(code, X, Yexp, model,opt.pg, opt.emul,opt.sim)
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
                                biais <- rnorm(1,0,sqrt(Cov))
                              } else
                              {
                                biais <- rnorm(n=self$n,0,sqrt(Cov))
                              }
                            } else
                            {
                              biais <- mvrnorm(n=self$n,rep(0,long),Cov)
                              biais <- apply(biais,1,mean)
                            }
                            return(list(biais=biais,cov=Cov))
                          },
                          pred = function(theta,thetaD,var,x.new)
                          {
                            self$disc <- self$discrepancy(theta,thetaD,var,x.new)
                            foo <- self$predTemp(theta,var,x.new)
                            y <- foo$y
                            yc  <- foo$yc
                            return(list(y=self$disc$biais+y,cov=self$disc$cov,yc=yc))
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
                            return(list(y=self$disc$biais+y,Cov.D=self$disc$cov,yc=yc,
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
                   # cat(paste("Mean of the biais:",round(mean(self$disc$biais),5),"\n",sep=" "))
                   # cat(paste("Covariance of the biais:",round(mean(self$disc$cov),5),"\n",sep=" "))
                   cat("Chosen kernel:", self$opt.disc$kernel.type)
                 }
)

