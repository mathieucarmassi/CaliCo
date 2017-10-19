#' A Reference Class to generates differents model objects
#'
#' @description See the function blabla which produces an instance of this class
#' This class comes with a set of methods, some of them being useful for the user:
#' See the documentation for blabla... Other methods
#'  should not be called as they are designed to be used during the optimization process.
#'
#' Fields should not be changed or manipulated by the user as they are updated internally
#' during the estimation process.
#' @field code a function which takes in entry X and theta
#' @field X the matrix of the forced variables
#' @field theta the vector of parameters to estimate
#' @field Yexp the experimental output
#' @field n the number of experiments
#' @field model the model choice see documentation
model.class <- R6::R6Class(classname = "model.class",
                 public = list(
                   code      = NULL,
                   X         = NULL,
                   Yexp      = NULL,
                   n         = NULL,
                   d         = NULL,
                   binf      = NULL,
                   bsup      = NULL,
                   emul.list = NULL,
                   model     = NULL,
                   initialize = function(code=NA,X=NA,Yexp=NA,model=NA,
                                        binf=NA,bsup=NA, emul.list=list(p=NA,n.emul=NA,PCA=NA))
                   {
                     self$code  <- code
                     self$X     <- X
                     self$Yexp  <- Yexp
                     self$n     <- length(Yexp)
                     self$binf  <- binf
                     self$bsup  <- bsup
                     if (is.null(dim(X)))
                     {
                       self$d       <- 1
                     } else
                     {
                       self$d       <- dim(X)[2]
                     }
                     self$emul.list <- emul.list
                     self$model     <- model
                     private$checkModels()
                     private$checkEmul()
                     #private$checkCode()
                     private$loadPackages()
                   }
                 ))

model.class$set("private","checkModels",
        function()
        {
          if (self$model != "model1" & self$model != "model2" & self$model != "model3" & self$model != "model4")
          {
            stop('Please elect a correct model')
          }
        })

model.class$set("private","loadPackages",
                function()
                {
                  library(R6)
                  library(DiceDesign)
                  library(DiceKriging)
                  library(FactoMineR)
                  library(Rcpp)
                  library(RcppArmadillo)
                  library(MASS)
                })

model.class$set("private","checkEmul",
                function()
                  {
                  N <- c("p","n.emul","PCA")
                  N2 <- names(self$emul.list)
                  for (i in 1:length(N))
                  {
                    if(names(self$emul.list)[i] != N[i])
                    {
                      stop(paste(N[i],"value is missing, please enter a correct value",sep=" "))
                    }
                  }
                })

model.class$set("private","checkCode",
                function()
                {
                  if (is.na(self$code))
                  {stop("Please enter a valid code")}
                })


model1.class <- R6::R6Class(classname = "model1.class",
                        inherit = model.class,
                        public=list(
                        m.exp = NULL,
                        V.exp = NULL,
                        initialize=function(code=NA, X=NA, Yexp=NA, model=NA, binf=NA, bsup=NA)
                        {
                          super$initialize(code, X, Yexp, model, binf, bsup)
                        },
                        fun = function(theta,sig2,Newdata)
                        {
                          y  <- self$code(Newdata,theta)+rnorm(self$n,0,sqrt(sig2))
                          yc <- self$code(Newdata,theta)
                          return(list(y=y, yc=yc))
                        },
                        likelihood = function(theta,sig2)
                        {
                          self$m.exp = self$code(self$X,theta)
                          self$V.exp = sig2*diag(self$n)
                          return(1/((2*pi)^(self$n/2)*det(self$V.exp)^(1/2))*exp(-1/2*t(self$Yexp-self$m.exp)%*%
                                                          invMat(self$V.exp)%*%(self$Yexp-self$m.exp)))
                        }
                        )
                        )


model1.class$set("public","plot",
                 function(theta,sig,Newdata)
                 {
                   res <- self$fun(theta,sig,Newdata)
                   gg.data <- data.frame(y=res$yc,x=seq(0,1,length.out=length(res$yc)),type="code result")
                   gg.data.noisy <- data.frame(y=res$y,x=seq(0,1,length.out=length(res$yc)),
                                             type="noisy")
                   gg.data.exp  <- data.frame(y=self$Yexp,x=seq(0,1,length.out=length(res$yc)),
                                              type="experiments")
                   gg.data <- rbind(gg.data,gg.data.noisy,gg.data.exp)
                   p <- ggplot(gg.data)+
                     geom_line(aes(y=y,x=x,col=type))+
                     theme_light()+
                     ylab("")+xlab("")+
                     theme(legend.position=c(0.65,0.86),
                           legend.text=element_text(size = '15'),
                           legend.title=element_blank(),
                           legend.key=element_rect(colour=NA),
                           axis.text=element_text(size=20))
                   return(p)
                 })



model3.class <- R6::R6Class(classname = "model3.class",
                        inherit = model1.class,
                        public=list(
                          funTemp  = NULL,
                          temp     = NULL,
                          initialize=function(code=NA, X=NA, Yexp=NA, model=NA, binf=NA, bsup=NA)
                          {
                            super$initialize(code, X, Yexp, model,binf, bsup)
                            self$funTemp <- super$fun
                          },
                          discrepancy = function(theta,thetaD,sig2)
                          {
                            y   <- self$funTemp(theta,sig2,self$X)$y
                            z   <- self$Yexp - y
                            Cov <- matrix(nr=self$n,nc=self$n)
                            if (is.null(ncol(self$X)))
                            {
                              Cov <- matrix(nr=self$n,nc=self$n)
                              for (j in 1:self$n)
                              {
                                for (i in 1:self$n)
                                {
                                  Cov[i,j] <- thetaD[1]*exp(-1/2*(sum((self$X[i]-self$X[j])^2)/thetaD[2])^2)
                                }
                              }
                            } else
                            {
                              Cov <- matrix(nr=self$n,nc=self$n)
                              for (j in 1:self$n)
                              {
                                for (i in 1:self$n)
                                {
                                  Cov[i,j] <- thetaD[1]*exp(-1/2*(sum((self$X[i,]-self$X[j,])^2)/thetaD[2])^2)
                                }
                              }
                            }
                            p <- eigen(Cov)$vectors
                            e <- eigen(Cov)$values
                            if (all(e>0)){} else
                            {
                              e[which(e<0)] <- 1e-4
                            }
                            d <- diag(e)
                            Cov <- t(p)%*%d%*%p
                            biais <- mvrnorm(n=100,rep(0,length(self$X)),Cov)
                            biais <- apply(biais,1,mean)
                            return(list(biais=biais,cov=Cov))
                          },
                          fun = function(theta,thetaD,sig2,Newdata)
                          {
                            res <- self$discrepancy(theta,thetaD,sig2)
                            foo <- self$funTemp(theta,sig2,Newdata)
                            y   <- foo$y
                            yc  <- foo$yc
                            return(list(y=res$biais+y,cov=res$cov,yc=yc))
                          }
                          )
)


model3.class$set("public","plot",
                 function(theta,thetaD,sig,Newdata)
                 {
                   res <- self$fun(theta,thetaD,sig,Newdata)
                   gg.data <- data.frame(y=res$yc,x=seq(0,1,length.out=length(res$yc)),type="code result")
                   gg.data.noisy <- data.frame(y=res$y,x=seq(0,1,length.out=length(res$yc)),
                                               type="with discrepancy and noise")
                   gg.data.exp  <- data.frame(y=self$Yexp,x=seq(0,1,length.out=length(res$yc)),
                                              type="experiments")
                   gg.data <- rbind(gg.data,gg.data.noisy,gg.data.exp)
                   p <- ggplot(gg.data)+
                     geom_line(aes(y=y,x=x,col=type))+
                     theme_light()+
                     ylab("")+xlab("")+
                     theme(legend.position=c(0.65,0.86),
                           legend.text=element_text(size = '15'),
                           legend.title=element_blank(),
                           legend.key=element_rect(colour=NA),
                           axis.text=element_text(size=20))
                   return(p)
                 })


model3.class$set("public","likelihood",
                 function(theta,thetaD,sig2,Newdata)
                 {
                   self$m.exp <- self$code(Newdata,theta)
                   temp <- self$fun(theta,thetaD,sig2,Newdata)
                   self$V.exp <- sig2*diag(self$n) + temp$cov
                   return(1/((2*pi)^(self$n/2)*det(self$V.exp)^(1/2))*exp(-1/2*t(self$Yexp-self$m.exp)%*%
                                                                invMat(self$V.exp)%*%(self$Yexp-self$m.exp)))
                 })



model2.class <- R6::R6Class(classname = "model2.class",
                        inherit = model.class,
                        public = list(
                          n.emul = NULL,
                          GP     = NULL,
                          DOE    = NULL,
                          p      = NULL,
                          binf   = NULL,
                          bsup   = NULL,
                          PCA    = NULL,
                          m.exp  = NULL,
                          V.exp  = NULL,
                        initialize = function(code=NA, X=NA, Yexp=NA, model=NA,opt.emul=NA,binf=NA,bsup=NA)
                        {
                          super$initialize(code, X, Yexp, model, binf, bsup, opt.emul)
                          self$binf   <- binf
                          self$bsup   <- bsup
                          self$n.emul <- opt.emul$n.emul
                          self$p      <- opt.emul$p
                          self$PCA    <- opt.emul$PCA
                          foo         <- self$surrogate()
                          self$GP     <- foo$GP
                          self$DOE    <- foo$DOE
                          print("The surrogate has been set up, you can now use the function")
                        },
                        surrogate = function()
                        {
                          Xcr <- scale(X)
                          V   <- attr(Xcr,"scaled:scale")
                          M   <- attr(Xcr,"scaled:center")
                          Dim <- self$p+self$d
                          if (self$PCA==TRUE)
                          {
                            D <- self$PCA.fun(X=self$X,Dim=Dim,n=self$n.emul,p=self$p,d=self$d,
                                              binf=self$binf,bsup=self$bsup,M=M,V=V)
                          } else
                          {
                            doe <- lhsDesign(self$n.emul,Dim)$design
                            doe <- maximinSA_LHS(doe)
                            doe <- doe$design
                            ### Getting back the value of the parameter generated by the DOE
                            binf.X <- apply(Xcr,2,min)
                            bsup.X <- apply(Xcr,2,max)
                            DOE <- unscale(doe[,1:self$d],binf.X,bsup.X)
                            if (self$d==1)
                            {
                              for (i in 1:self$n.emul)
                              {
                                if (self$d==1)
                                {
                                  DOE[i] <- DOE[i]*V+M
                                } else
                                {
                                  DOE[i] <- DOE[i]*V[i]+rep(M[i],self$n.emul)
                                }
                              }
                            } else {
                              for (i in 1:self$n.emul)
                              {
                                DOE[,i] <- DOE[,i]*V[i]+rep(M[i],self$n.emul)
                              }
                            }
                            doeParam <- unscale(doe[,(Dim-self$p+1):Dim],self$binf,self$bsup)
                            ### Matrix D contains the final value for the DOE
                            D <- cbind(DOE,doeParam)
                          }
                          ### Generating the response
                          z <- self$code(D[,1:(Dim-self$p)],D[,(Dim-self$p+1):Dim])
                          ### Converting D as a data.frame for the km function
                          D <- as.data.frame(D)
                          #colnames(D) <- c("V1","V2","V3","V4","V5","V6")
                          ### Creation of the Gaussian Process with estimation of hyperpameters
                          GP <- km(formula =~1, design=D, response = z,covtype = "matern5_2")
                          return(list(GP=GP,DOE=D))
                        },
                        fun = function(theta,sig2,Newdata)
                        {
                          # options(warn=-1)
                          if(self$p==1)
                          {
                            Xnew <- cbind(Newdata,rep(theta,self$n))
                          } else
                          {
                            Xtemp <- matrix(rep(theta,c(self$n,self$n)),nr=self$n,nc=self$p)
                            Xnew  <- cbind(Newdata,Xtemp)
                          }
                          Xnew <- as.data.frame(Xnew)
                          names(Xnew) <- c("DOE","doeParam")
                          pr <- predict(self$GP,newdata=as.data.frame(Xnew),type="UK",
                                        cov.compute=TRUE,interval="confidence")
                          err <- rnorm(n=self$n,mean = 0,sd=sqrt(sig2))
                          return(list(y=pr$mean+err,Cov.GP=pr$cov,yc=pr$mean,lower=pr$lower95,upper=pr$upper95))
                        })
                        )

model2.class$set("public","PCA.fun",
                function(X,Dim,n,p,d,binf,bsup,M,V)
                {
                  PCA.sim <- PCA(X,graph = FALSE)
                  ### Coordinates in the new uncorrelated space of the initial points
                  B <- PCA.sim$ind$coord
                  ### Transition matrix
                  P <- sqrt(PCA.sim$var$contrib[,1:d])/10 * sign(PCA.sim$var$coord[,1:d])
                  ### Establishment of the DOE
                  doe <- lhsDesign(n,Dim)$design
                  doe <- maximinSA_LHS(doe)
                  doe <- doe$design
                  ### Boundaries of the points in the new coordinates
                  binf.X <- apply(B,2,min)
                  bsup.X <- apply(B,2,max)
                  ### Unscaling the doe into the right bounds
                  DOE <- unscale(doe[,1:d],binf.X,bsup.X)
                  ### Matrix containing the points from the DOE but in the initial coordinates
                  A <- t(P%*%t(DOE))
                  ### Multiplication of each components by the variance and add the mean
                  for (i in 1:d)
                  {
                    A[,i] <- A[,i]*V[i]+rep(M[i],n)
                  }
                  ### Getting back the value of the parameter generated by the DOE
                  doeParam <- unscale(doe[,(Dim-p+1):Dim],binf,bsup)
                  ### Matrix D contains the final value for the DOE
                  D <- cbind(A,doeParam)
                  return(D)
                })

model2.class$set("public","likelihood",
                 function(theta,sig2,Newdata)
                 {
                   temp <- self$fun(theta,sig2,Newdata)
                   self$m.exp <- temp$yc
                   self$V.exp <- sig2*diag(self$n) + temp$Cov.GP
                   return(1/((2*pi)^(self$n/2)*det(self$V.exp)^(1/2))*exp(-1/2*t(self$Yexp-self$m.exp)%*%
                                                                  invMat(self$V.exp)%*%(self$Yexp-self$m.exp)))
                 })


model2.class$set("public","plot",
                 function(theta,sig,Newdata,points=TRUE)
                 {
                   res <- self$fun(theta,sig,Newdata)
                   gg.data <- data.frame(y=res$yc,x=seq(0,1,length.out=length(res$yc)),
                                          lower=res$lower,upper=res$upper,type="Gaussian Process",
                                         fill="90% credibility interval")
                   gg.data.exp <- data.frame(y=self$Yexp,x=seq(0,1,length.out=length(res$yc)),lower=res$lower,
                                            upper=res$upper,type="experiment",fill="90% credibility interval")
                   gg.data <- rbind(gg.data,gg.data.exp)
                   gg.points <- data.frame(x=self$DOE[,1],y=self$code(self$DOE[,1],theta))
                   p <- ggplot(gg.data)+ geom_ribbon(aes(ymin=lower,ymax=upper,x=x,fill=fill),alpha=0.3)+
                     geom_line(aes(y=y,x=x,col=type))+
                     theme_light()+
                     ylab("")+xlab("")+
                     scale_fill_manual("",values=c("grey12"))+
                     theme(legend.position=c(0.65,0.86),
                           legend.text=element_text(size = '15'),
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
                 })



model4.class <- R6::R6Class(classname = "model4.class",
                        inherit = model2.class,
                        public=list(
                          funC = NULL,
                          initialize=function(code=NA, X=NA, Yexp=NA, model=NA,opt.emul=NA,binf=NA,bsup=NA)
                          {
                            super$initialize(code, X, Yexp, model, opt.emul, binf, bsup)
                            self$funC <- super$fun
                          },
                          discrepancy = function(theta,thetaD,sig2)
                          {
                            y   <- self$funTemp(theta,sig2,self$X)$y
                            z   <- self$Yexp - y
                            Cov <- matrix(nr=self$n,nc=self$n)
                            if (is.null(ncol(self$X)))
                            {
                              Cov <- matrix(nr=self$n,nc=self$n)
                              for (j in 1:self$n)
                              {
                                for (i in 1:self$n)
                                {
                                  Cov[i,j] <- thetaD[1]*exp(-1/2*(sum((self$X[i]-self$X[j])^2)/thetaD[2])^2)
                                }
                              }
                            } else
                            {
                              Cov <- matrix(nr=self$n,nc=self$n)
                              for (j in 1:self$n)
                              {
                                for (i in 1:self$n)
                                {
                                  Cov[i,j] <- thetaD[1]*exp(-1/2*(sum((self$X[,i]-self$X[,j])^2)/thetaD[2])^2)
                                }
                              }
                            }
                            p <- eigen(Cov)$vectors
                            e <- eigen(Cov)$values
                            if (all(e>0)){} else
                            {
                              e[which(e<0)] <- 1e-4
                            }
                            d <- diag(e)
                            Cov <- t(p)%*%d%*%p
                            biais <- mvrnorm(n=100,rep(0,length(self$X)),Cov)
                            biais <- apply(biais,1,mean)
                            return(list(biais=biais,cov=Cov))
                          },
                          fun = function(theta,thetaD,sig2,Newdata)
                          {
                            print(self$discrepancy)
                            res <- self$discrepancy(theta,thetaD,sig2)
                            foo <- self$funC(theta,sig2,Newdata)
                            y <- foo$y
                            Cov.GP <- foo$cov
                            yc <- foo$yc
                            lower <- foo$lower
                            upper <- foo$upper
                            return(list(y=res$biais+y,Cov.D=res$cov,yc=yc,lower=lower,upper=upper,Cov.GP=Cov.GP))
                          })
)


model4.class$set("public","likelihood",
                 function(theta,thetaD,sig2,Newdata)
                 {
                   temp <- self$fun(theta,thetaD,sig2,Newdataœ)
                   self$m.exp <- temp$yc
                   self$V.exp <- sig2*diag(self$n) + temp$Cov.GP +temp$Cov.D
                   return(1/((2*pi)^(self$n/2)*det(self$V.exp)^(1/2))*exp(-1/2*t(self$Yexp-self$m.exp)%*%
                                                                  invMat(self$V.exp)%*%(self$Yexp-self$m.exp)))
                 })



model4.class$set("public","plot",
                 function(theta,sig,Newdata,points=TRUE)
                 {
                   res <- self$fun(theta,thetaD,sig,Newdata)
                   gg.data <- data.frame(y=res$yc,x=seq(0,1,length.out=length(res$yc)),
                                         lower=res$lower,upper=res$upper,type="Gaussian Process",
                                         fill="90% credibility interval")
                   gg.data.dis <- data.frame(y=res$y,x=seq(0,1,length.out=length(res$yc)),lower=res$lower,
                                             upper=res$upper,type="experiment",fill="90% credibility interval")
                   gg.data.exp <- data.frame(y=self$Yexp,x=seq(0,1,length.out=length(res$yc)),lower=res$lower,
                                             upper=res$upper,type="experiment",fill="90% credibility interval")
                   gg.data <- rbind(gg.data,gg.data.dis,gg.data.exp)
                   gg.points <- data.frame(x=self$DOE[,1],y=self$code(self$DOE[,1],theta))
                   p <- ggplot(gg.data)+ geom_ribbon(aes(ymin=lower,ymax=upper,x=x,fill=fill),alpha=0.3)+
                     geom_line(aes(y=y,x=x,col=type))+
                     theme_light()+
                     ylab("")+xlab("")+
                     scale_fill_manual("",values=c("grey12"))+
                     theme(legend.position=c(0.65,0.86),
                           legend.text=element_text(size = '15'),
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
                 })

