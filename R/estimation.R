#' A Reference Class to generate different calibration methods after generating a model from
#' model.class
#'
#' @examples
#'
#' @export
estim.class <- R6::R6Class(classname = "estim.class",
                 public = list(
                    code        = NULL,
                    X           = NULL,
                    Yr          = NULL,
                    Yexp        = NULL,
                    model       = NULL,
                    type.prior  = NULL,
                    opt.emul    = NULL,
                    opt.prior   = NULL,
                    opt.estim   = NULL,
                    burnIn      = NULL,
                    binf        = NULL,
                    bsup        = NULL,
                    logTest.fun = NULL,
                    md          = NULL,
                    mdCV        = NULL,
                    pr          = NULL,
                    out         = NULL,
                    outCV       = NULL,
                    type.valid  = NULL,
                    opt.valid   = NULL,
                    initialize = function(code=NA,X=NA,Yexp=NA,model=NA,type.prior=NA,
                                          opt.emul=NA,opt.prior=NA,opt.estim=NA,type.valid=NA,opt.valid=NA)
                    {
                      self$code          <- code
                      self$X             <- X
                      self$Yexp          <- Yexp
                      self$model         <- model
                      self$opt.emul      <- opt.emul
                      self$opt.prior     <- opt.prior
                      self$opt.estim     <- opt.estim
                      self$burnIn        <- 0.1*opt.estim$Nmh
                      self$type.prior    <- type.prior
                      self$type.valid    <- type.valid
                      self$opt.valid     <- opt.valid
                      private$checkValid()
                      private$checkup()
                      self$pr            <- prior(self$type.prior,opt.prior,log=TRUE)
                      self$binf          <- private$boundaries()$binf
                      self$bsup          <- private$boundaries()$bsup
                      bound <- private$boundaries()
                      bound <- list(binf=bound[[1]][1:length(bound[[1]])-1],bsup=bound[[2]][1:length(bound[[2]])-1])
                      self$opt.emul      <- c(self$opt.emul,bound,list(DOE=NULL))
                      self$Yr            <- self$code(self$X,self$opt.estim$thetaInit)
                      if (is.null(self$type.valid))
                      {
                        self$md <- model(code,X,Yexp,model,self$opt.emul)
                        self$logTest.fun   <- self$logLikelihood(model)
                        self$out           <- self$estimation(self$md,self$Yexp)
                      } else
                      {
                        cat("#############################################\n")
                        cat("##### --- Begin of the calibration --- ######\n")
                        cat("#############################################\n")
                        self$md <- model(code,X,Yexp,model,opt.emul,binf=self$binf[1],
                                         bsup=self$bsup[1])
                        self$logTest.fun   <- self$logLikelihood(model)
                        self$out           <- self$estimation(self$md,self$Yexp)
                        cat("###################################################\n")
                        cat("##### --- End of the regular calibration --- ###### \n")
                        cat("###################################################\n")
                        self$outCV <- self$validation(self$type.valid,self$opt.valid)
                      }
                    },
                    estimation = function(md,Yexp)
                    {
                      MCMC <- function(model)
                      {
                        switch(model,
                               model1={return(MetropolisHastingsCpp)},
                               model2={return(MetropolisHastingsCpp)},
                               model3={return(MetropolisHastingsCppD)},
                               model4={return(MetropolisHastingsCppD)}
                        )
                      }
                      MetropolisCpp <- MCMC(self$model)
                      out <- MetropolisCpp(self$md$fun,self$opt.estim$Ngibbs,
                                                   self$opt.estim$Nmh,self$opt.estim$thetaInit,
                                                   self$opt.estim$k,self$opt.estim$sig,Yexp,
                                                   self$binf,self$bsup,self$logTest.fun)
                      return(out)
                    }
                   ))


estim.class$set("public","logLikelihood",
  function(model)
  {
    switch(model,
           model1={return(self$logTest)},
           model2={return(self$logTest)},
           model3={return(self$logTestD)},
           model4={return(self$logTestD)}
    )
  }
)


estim.class$set("private","boundaries",
                function()
                {
                  binf <- self$pr[[1]]$binf
                  bsup <- self$pr[[1]]$bsup
                  for (i in 2:length(self$type.prior))
                  {
                    binf <- c(binf,self$pr[[i]]$binf)
                    bsup <- c(bsup,self$pr[[i]]$bsup)
                  }
                  return(list(binf=binf,bsup=bsup))
                })


estim.class$set("private","checkValid",
                function()
                {
                  if (is.null(self$type.valid)==FALSE)
                  {
                    if(is.null(self$opt.valid)==FALSE)
                    {
                      if (self$type.valid != "loo" & self$type.valid != "kfold")
                      {
                        stop("Plese select a validation method ('loo' or 'kfold')")
                      }
                      if (self$type.valid == "loo" & is.null(self$opt.valid$n.CV))
                      {
                        stop("You have selected the leave one out validation, please enter a number of repetition in opt.valid")
                      }
                      if (self$type.valid == "kfold" & is.null(self$opt.valid$n.CV))
                      {
                        stop("You have selected the k-fold validation, please enter a number of repetition in opt.valid")
                      }
                      if (self$type.valid == "kfold" & is.null(self$opt.valid$k))
                      {
                        stop("You have selected the k-fold validation, please enter a valid k in opt.valid")
                      }
                    }
                  }
                })


estim.class$set("public","logTest",
                function(theta,sig2)
                {
                  if (length(self$type.prior) == 1)
                  {
                      return(log(self$md$likelihood(theta,sig2))+self$pr$prior(theta))
                  } else
                  {
                      s <- 0
                      for (i in 1:(length(theta)))
                      {
                        s <- s + self$pr[[i]]$prior(theta[i])
                      }
                      s <- s + self$pr[[(length(theta)+1)]]$prior(sig2)
                      return(log(self$md$likelihood(theta,sig2)) + s)
                  }
                })

estim.class$set("public","logTestD",
                function(theta,thetaD,sig2)
                {
                  s <- 0
                  for (i in 1:(length(theta)))
                  {
                    s <- s + self$pr[[i]]$prior(theta[i])
                  }
                  for (j in 1:(length(thetaD)))
                  {
                    s <- s + self$pr[[length(theta)+j]]$prior(thetaD[j])
                  }
                  s <- s + self$pr[[(length(theta)+1)]]$prior(sig2)
                  return(log(self$md$likelihood(theta,thetaD,sig2)) + s)
                })


estim.class$set("private","checkup",
                function()
                {
                 N <- c("Ngibbs","Nmh","thetaInit","k","sig")
                 for (i in 1:length(N))
                 {
                   if(names(self$opt.estim)[i]!=N[i])
                   {
                     stop(paste(N[i],"value is missing in opt.estim",sep=" "))
                   }
                   if(is.na(names(self$opt.estim)[i]))
                   {
                     stop(paste(N[i],"value is not correct, please enter a right value",sep=" "))
                   }
                 }}
		)

estim.class$set("public","plot",
               function(graph=c("acf","chains","densities","output"),separated=FALSE,CI=TRUE,select.X=NULL)
                 {
                    n      <- length(self$type.prior)
                    a      <- list()
                    m      <- list()
                    p      <- list()
                    output <- function()
                    {
                      if (self$md$model=="model1" | self$md$model=="model2")
                      {
                        return(self$plotComp(CI,select.X))
                      }else
                      {
                        return(self$plotCompD(CI,select.X))
                      }
                    }
                    for (i in 1:n)
                    {
                      dplot2  <- data.frame(data=self$out$THETA[-c(1:self$burnIn),i],type="posterior")
                      a[[i]]  <- self$acf(i)
                      m[[i]]  <- self$mcmc(i)
                      p[[i]]  <- self$pr[[i]]$plot() + geom_density(data=dplot2,kernel="gaussian",adjust=3,alpha=0.1)
                    }
                    if (separated==TRUE)
                    {
                      if (length(graph)==1)
                      {
                        if (graph=="acf"){return(a)}
                        if (graph=="chains"){return(m)}
                        if (graph=="densities"){return(p)}
                        if (graph=="output"){return(output())}
                      }
                      if (length(graph)>1)
                      {
                        if (graph==c("acf","chains")){a;m}
                        if (graph==c("acf","densities")){a;p}
                        if (graph==c("acf","output")){a;output()}
                        if (graph==c("chains","densities")){m;p}
                        if (graph==c("chains","output")){m;output()}
                        if (graph==c("densities","output")){p;output()}
                      }
                      if (length(graph)>2)
                      {
                        if (graph==c("acf","chains","densities")){a;m;p}
                        if (graph==c("acf","chains","output")){a;m;output()}
                        if (graph==c("acf","densities","output")){a;p;output()}
                        if (graph==c("chains","densities","output")){m;p;output()}
                      }
                      if (length(graph)>3)
                      {
                        if (graph==c("acf","chains","densities","output")){a;m;p;output}
                      }
                    } else{
                      if (length(graph)==1)
                      {
                        if (graph=="acf"){return(do.call(grid.arrange,a))}
                        if (graph=="chains"){return(do.call(grid.arrange,m))}
                        if (graph=="densities"){return(do.call(grid.arrange,p))}
                        if (graph=="output"){return(output())}
                      }
                      if (length(graph)>1)
                      {
                        if (graph==c("acf","chains")){do.call(grid.arrange,a);do.call(grid.arrange,m)}
                        if (graph==c("acf","densities")){do.call(grid.arrange,a);do.call(grid.arrange,p)}
                        if (graph==c("acf","output")){do.call(grid.arrange,a);output()}
                        if (graph==c("chains","densities")){do.call(grid.arrange,m);do.call(grid.arrange,p)}
                        if (graph==c("chains","output")){do.call(grid.arrange,m);output()}
                        if (graph==c("densities","output")){do.call(grid.arrange,p);output()}
                      }
                      if (length(graph)>2)
                      {
                        if (graph==c("acf","chains","densities")){do.call(grid.arrange,a);do.call(grid.arrange,m);
                          do.call(grid.arrange,p)}
                        if (graph==c("acf","chains","output")){do.call(grid.arrange,a);do.call(grid.arrange,m);output()}
                        if (graph==c("acf","densities","output")){do.call(grid.arrange,a);do.call(grid.arrange,p);output()}
                        if (graph==c("chains","densities","output")){do.call(grid.arrange,m);do.call(grid.arrange,p);output()}
                      }
                      if (length(graph)>3)
                      {
                        if (graph==c("acf","chains","densities","output")){do.call(grid.arrange,a);
                          do.call(grid.arrange,m);do.call(grid.arrange,p);output()}
                      }
                    }
               })

estim.class$set("public","gg",
                function(i)
                {
                  dplot  <- data.frame(data=self$out$THETA[-c(1:self$burnIn),i],type="posterior")
                  p <- ggplot(dplot,aes(data,fill=type,color=type)) +
                    geom_density(kernel = "gaussian",adjust=3,alpha=0.1)+
                    theme_light()+xlab("")+ylab("")+ xlim(self$binf[i],self$bsup[i])+
                    scale_fill_manual( values = "blue")+
                    scale_color_manual(values = "blue")+
                    theme(legend.position=c(0.86,0.86),
                          legend.text=element_text(face="bold",size = '20'),
                          legend.title=element_blank(),
                          legend.key=element_rect(colour=NA),
                          axis.text=element_text(size=20))+
                    geom_hline(aes(yintercept = 0))
                  return(p)
                })


estim.class$set("public","acf",
                function(i)
                  {
                    bacf   <- acf(self$out$THETA[-c(1:self$burnIn),i], plot = FALSE)
                    bacfdf <- with(bacf, data.frame(lag, acf))
                    p      <- ggplot(data = bacfdf, mapping = aes(x = lag, y = acf))+
                        geom_hline(aes(yintercept = 0))+
                        geom_segment(mapping = aes(xend = lag, yend = 0))+
                        xlab("")+ylab("")+theme_light()
                    return(p)
                  }
)


estim.class$set("public","mcmc",
                function(i)
                {
                  n <- length(self$out$THETA[-c(1:self$burnIn),i])
                  resgg <- data.frame(inc=c(1:n),data=self$out$THETA[-c(1:self$burnIn),i])
                  p   <- ggplot(data=resgg, aes(x=inc,y=data))+geom_line()+ylab("")+
                      xlab("")+theme_light()
                  return(p)
                }
)

estim.class$set("public","plotComp",
                function(CI=TRUE,depend.X=TRUE)
                {
                  m <- self$out$THETA[-c(1:self$burnIn),]
                  Dist <- matrix(nr=nrow(m),nc=length(self$Yexp))
                  Dim <- length(self$type.prior)
                  for (i in 1:nrow(m))
                  {
                    Dist[i,] <- self$md$fun(m[i,1:(Dim-1)],m[i,Dim])$y
                  }
                  lowerPost <- apply(Dist,2,quantile,probs=0.05)
                  upperPost <- apply(Dist,2,quantile,probs=0.95)
                  meanPost  <- apply(Dist,2,quantile,probs=0.5)
                  if (depend.X==FALSE)
                  {
                    X <- self$X[,1]
                  }else
                  {
                    if (is.null(dim(self$X)))
                    {
                      X <- self$X
                    } else{
                      if (dim(self$X)[2]>1)
                      {
                        stop("You an X dimension over 2, please desactivate the option depend.X to visualize your result")
                      }
                    }
                  }
                  if (CI==TRUE)
                  {
                    dplot <- data.frame(Y=meanPost,x=X,type='calibrated',lower=lowerPost,upper=upperPost,
                                        fill="90% credibility interval a posteriori")
                    dplot2 <- data.frame(Y=self$Yr, x=X, type='experiment',lower=lowerPost,upper=upperPost,
                                         fill="90% credibility interval a posteriori")
                    dplot <- rbind(dplot,dplot2)
                    p <- ggplot(dplot) + geom_line(aes(x=x,y=Y,color=type))+
                      geom_ribbon(aes(ymin=lower, ymax=upper, x=x,fill=fill), alpha = 0.3) +
                      scale_fill_manual("",values=c("blue4","blue4")) +
                      theme_light() + xlab("") + ylab("") +
                      theme(legend.position=c(0.65,0.86),
                            legend.text=element_text(size = '15'),
                            legend.title=element_blank(),
                            legend.key=element_rect(colour=NA),
                            axis.text=element_text(size=20))
                  } else
                  {
                    dplot  <- data.frame(x=X,data=meanPost,type="calibrated")
                    dplot2 <- data.frame(x=X,data=self$Yr,type="experiments")
                    dplot  <- rbind(dplot,dplot2)
                    p <- ggplot(data=dplot, aes(x=x,y=data,color=type))+geom_line()+ylab("")+xlab("")+
                      theme_light()+
                      theme(legend.position=c(0.86,0.86),
                            legend.text=element_text(face="bold",size = '20'),
                            legend.title=element_blank(),
                            legend.key=element_rect(colour=NA),
                            axis.text=element_text(size=20))
                  }
                    return(p)
                })


estim.class$set("public","plotCompD",
                function(CI=TRUE,depend.X=TRUE)
                {
                  m <- self$out$THETA[-c(1:self$burnIn),]
                  Dist <- matrix(nr=nrow(m),nc=length(self$Yexp))
                  Dim <- length(self$type.prior)
                  for (i in 1:nrow(m))
                  {
                    Dist[i,] <- self$md$fun(m[i,1:(Dim-3)],m[i,(Dim-2):(Dim-1)],m[i,Dim])$y
                  }
                  lowerPost <- apply(Dist,2,quantile,probs=0.05)
                  upperPost <- apply(Dist,2,quantile,probs=0.95)
                  meanPost  <- apply(Dist,2,quantile,probs=0.5)
                  if (depend.X==FALSE)
                  {
                    X <- self$X[,1]
                  }else
                  {
                    if (is.null(dim(self$X)))
                    {
                      X <- self$X
                    } else{
                      if (dim(self$X)[2]>1)
                      {
                        stop("You an X dimension over 2, please desactivate the option depend.X to visualize your result")
                      }
                    }
                  }
                  if (CI==TRUE)
                  {
                    dplot <- data.frame(Y=meanPost,x=self$X,type='calibrated',lower=lowerPost,upper=upperPost,
                                        fill="90% credibility interval a posteriori")
                    dplot2 <- data.frame(Y=self$Yr, x=self$X, type='experiment',lower=lowerPost,upper=upperPost,
                                         fill="90% credibility interval a posteriori")
                    dplot <- rbind(dplot,dplot2)
                    p <- ggplot(dplot) + geom_line(aes(x=x,y=Y,color=type))+
                      geom_ribbon(aes(ymin=lower, ymax=upper, x=x,fill=fill), alpha = 0.3) +
                      scale_fill_manual("",values=c("blue4","blue4")) +
                      theme_light() + xlab("") + ylab("") +
                      theme(legend.position=c(0.65,0.86),
                            legend.text=element_text(size = '15'),
                            legend.title=element_blank(),
                            legend.key=element_rect(colour=NA),
                            axis.text=element_text(size=20))
                  } else
                  {
                    dplot  <- data.frame(x=self$X,data=meanPost,type="calibrated")
                    dplot2 <- data.frame(x=self$X,data=self$Yr,type="experiments")
                    dplot  <- rbind(dplot,dplot2)
                    p <- ggplot(data=dplot, aes(x=x,y=data,color=type))+geom_line()+ylab("")+xlab("")+
                      theme_light()+
                      theme(legend.position=c(0.86,0.86),
                            legend.text=element_text(face="bold",size = '20'),
                            legend.title=element_blank(),
                            legend.key=element_rect(colour=NA),
                            axis.text=element_text(size=20))
                  }
                  return(p)
                })

estim.class$set("public","print",
                function()
                {
                  cat('the main results of the calibration are')
                }
)

estim.class$set("public","validation",
                function(type.valid,opt.valid)
                {
                  cat("\n")
                  cat("####################################################\n")
                  cat("##### --- Bigining of the cross validation --- #####\n")
                  cat("####################################################\n")
                  id       <- matrix(nr=opt.valid$n.CV,nc=1)
                  # Generates differents lists which store results from calibration (Rescal),
                  # the quantiles and the coverTau
                  ResCal   <- q5 <- m <- q95  <- list()
                  RMSE     <- coverTau  <- rep(0,opt.valid$n.CV)
                  qtheta5  <- qtheta95 <- mtheta <- list()
                  Yc       <- matrix(nr=opt.valid$n.CV,nc=1)
                  # Storage of the initial data because modification of self$X will be operated later on
                  dataInit <- self$X
                  # Storage of the initial output because modification of self$Yexp will be operated later on
                  YexpInit <- self$Yexp
                  # Beginning of the cross validation
                  for (i in 1:opt.valid$n.CV)
                  {
                    cat("####################################################\n")
                    cat(paste("##### --- iteration number --- #### ",i," #### --- #####\n",sep=""))
                    cat("####################################################\n")
                    if (type.valid=="loo")
                    {
                      slt <- sample(self$X,1)
                    } else
                    {
                      if (type.valid == "kfold")
                      {
                        slt <- sample(self$X,opt.valid$k)
                      }
                    }
                    # Getting back calibration and validation data set
                    id[i]             <- which(slt==dataInit)
                    dataCalib         <- dataInit[-id[i]]
                    dataValid         <- dataInit[id[i]]
                    Ycalib            <- YexpInit[-id[i]]
                    Yvalid            <- YexpInit[id[i]]
                    # Generate a model.class objet from calibration data set
                    self$mdCV         <- model(self$code,dataCalib,Ycalib,self$model,self$opt.emul,binf=self$binf[1],
                                               bsup=self$bsup[1])
                    # Generate a liklihood from the selected model
                    self$X            <- dataCalib
                    self$Yexp         <- Ycalib
                    self$logTest.fun  <- self$logLikelihood(self$model)
                    yres              <- matrix(nr=self$opt.estim$Nmh-self$burnIn,nc=length(Ycalib))
                    # MCMC for the iteration i. The results are stored in ResCal
                    ResCal[[i]]       <- self$estimation(self$mdCV,Ycalib)$THETA
                    FlushCPP()
                    # DimC represents the number of parameters to calibrate
                    DimC              <- ncol(ResCal[[i]])
                    # Construction of the posterior distribution for each points found in the MCMC
                    for (j in self$burnIn:(self$opt.estim$Nmh-1))
                    {
                      inc        <- j-self$burnIn+1
                      yres[inc,] <- self$mdCV$fun(ResCal[[i]][j,c(1:(DimC-1))],ResCal[[i]][j,DimC])$y
                    }
                    # To let the user the acces to the 90% confidence interval a posteriori
                    q        <- apply(yres,2,quantile,probs=c(0.05,0.5,0.95))
                    q5[[i]]  <- q[1,]
                    m[[i]]   <- q[2,]
                    q95[[i]] <- q[3,]
                    # Get the estimation in the posterior distribution
                    qtheta        <- apply(ResCal[[i]],2,quantile,probs=c(0.05,0.5,0.95))
                    qtheta5[[i]]  <- qtheta[1,]
                    mtheta[[i]]   <- qtheta[2,]
                    qtheta95[[i]] <- qtheta[3,]
                    # Generate a model.class objet from validation data set
                    self$mdCV         <- model(self$code,dataValid,Yvalid,self$model,self$opt.emul,
                                               binf=mean(qtheta5[[i]]),bsup=mean(qtheta95[[i]]))
                    # Generate a liklihood from the selected model
                    self$X            <- dataValid
                    self$Yexp         <- Yvalid
                    self$logTest.fun  <- self$logLikelihood(self$model)
                    # Compute the 90% confidence interval for the validation
                    yvalid            <- self$mdCV$fun(mtheta[[i]][c(1:(DimC-1))],mtheta[[i]][DimC])$y
                    qvalid5           <- self$mdCV$fun(qtheta5[[i]][c(1:(DimC-1))],qtheta5[[i]][DimC])$y
                    qvalid95          <- self$mdCV$fun(qtheta95[[i]][c(1:(DimC-1))],qtheta95[[i]][DimC])$y
                    # For each iteration of the CV check if the point is in the 90% confidence interval
                    # if loo is activated
                    if (Yvalid > qvalid5 & Yvalid < qvalid95)
                    {
                      coverTau[i] <- coverTau[i] + 1
                    }
                    RMSE[i] <- sqrt(mean((yvalid-Yvalid)^2))
                  }
                  self$X    <- dataInit
                  self$Yexp <- YexpInit
                  return(list(coverTau=mean(coverTau),RMSE=mean(RMSE)))
                }
)



