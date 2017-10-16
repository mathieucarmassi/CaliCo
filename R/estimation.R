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
                    log         = NULL,
                    opt.emul    = NULL,
                    opt.prior   = NULL,
                    opt.estim   = NULL,
                    burnIn      = NULL,
                    binf        = NULL,
                    bsup        = NULL,
                    logTest.fun = NULL,
                    md          = NULL,
                    pr          = NULL,
                    out         = NULL,
                    initialize = function(code=NA,X=NA,Yr=NA,Yexp=NA,model=NA,type.prior=NA,log=TRUE,
                                          opt.emul=NA,opt.prior=NA,opt.estim=NA)
                    {
                      self$code          <- code
                      self$X             <- X
                      self$Yr            <- Yr
                      self$Yexp          <- Yexp
                      self$model         <- model
                      self$log           <- log
                      self$opt.emul      <- opt.emul
                      self$opt.prior     <- opt.prior
                      self$opt.estim     <- opt.estim
                      self$burnIn        <- 0.1*opt.estim$Nmh
                      self$type.prior    <- type.prior
                      private$checkup()
                      self$pr            <- prior(self$type.prior,opt.prior,log=TRUE)
                      self$binf          <- private$boundaries()$binf
                      self$bsup          <- private$boundaries()$bsup
                      self$md            <- model(code,X,Yexp,model,opt.emul,binf=self$binf[1:(length(self$type.prior)-1)],
                                                  bsup=self$bsup[1:(length(self$type.prior)-1)])
                      self$logTest.fun   <- self$logTest
                      self$out           <- self$estimation()
                    },
                    estimation = function()
                    {
                      # print(self$md$fun)
                      # print(self$opt.estim$Ngibbs)
                      # print(self$opt.estim$Nmh)
                      # print(self$opt.estim$thetaInit)
                      # print(self$opt.estim$k)
                      # print(self$opt.estim$sig)
                      # print(self$Yexp)
                      # print(self$binf)
                      # print(self$bsup)
                      # print(self$logTest.fun)
                      out <- MetropolisHastingsCpp(self$md$fun,self$opt.estim$Ngibbs,
                                                   self$opt.estim$Nmh,self$opt.estim$thetaInit,
                                                   self$opt.estim$k,self$opt.estim$sig,self$Yexp,
                                                   self$binf,self$bsup,self$logTest.fun)
                      return(out)
                    }
                   ))

estim.class$set("private","boundaries",
                function()
                {
                  binf <- self$pr[[1]]$binf
                  bsup <- self$pr[[1]]$bsup
                  for (i in 2:length(self$type.prior))
                  {
                    binf <- cbind(binf,self$pr[[i]]$binf)
                    bsup <- cbind(bsup,self$pr[[i]]$bsup)
                  }
                  return(list(binf=binf,bsup=bsup))
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
                 }
                })


estim.class$set("public","plot",
               function(separated=FALSE,CI=TRUE)
                 {
                    n      <- length(self$type.prior)
                    a      <- list()
                    m      <- list()
                    p      <- list()
                    for (i in 1:n)
                    {
                      dplot2  <- data.frame(data=self$out$THETA[-c(1:self$burnIn),i],type="posterior")
                      a[[i]]  <- self$acf(i)
                      m[[i]]  <- self$mcmc(i)
                      p[[i]]  <- self$pr[[i]]$plot() + geom_density(data=dplot2,kernel="gaussian",adjust=3,alpha=0.1)
                    }
                    if (separated==TRUE)
                    {
                      a
                      m
                      p
                    } else{
                      do.call(grid.arrange,a)
                      do.call(grid.arrange,m)
                      do.call(grid.arrange,p)
                    }
                    self$plotComp(CI)
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
                function(CI=TRUE)
                  {
                    if (CI==TRUE)
                    {
                      m <- apply(self$out$THETA[-c(1:self$burnIn),],2,mean)
                      quantilesPost <- apply(self$out$THETA[-c(1:self$burnIn),],2,quantile,probs=c(0.05,0.95,0.5))
                      lowerPost <- self$code(self$X,quantilesPost[1,1:(length(m)-1)])
                      upperPost <- self$code(self$X,quantilesPost[2,1:(length(m)-1)])
                      dplot <- data.frame(Y=self$code(self$X,m[-length(m)]),x=self$X,type='calibrated',lower=lowerPost,upper=upperPost,
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
                      m      <- apply(self$out$THETA[-c(1:self$burnIn),],2,mean)
                      dplot  <- data.frame(x=self$X,data=self$code(m[-length(m)]),type="calibrated")
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


