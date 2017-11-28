#' A Reference Class to generate different calibration methods after generating a model from
#' model.class
#'
#' @examples
#'
#' @export
calibrate.class <- R6::R6Class(classname = "calibrate.class",
                           public = list(
                             md        = NULL,
                             pr        = NULL,
                             opt.estim = NULL,
                             opt.valid = NULL,
                             logPost   = NULL,
                             mcmc      = NULL,
                             output    = NULL,
                             errorCV   = NULL,
                             initialize = function(md=NA,pr=NA,opt.estim=NA,opt.valid=NULL)
                             {
                               library(parallel)
                               self$md        <- md
                               self$pr        <- pr
                               self$opt.estim <- opt.estim
                               self$opt.valid <- opt.valid
                               self$logPost   <- private$logLikelihood(self$md$model)
                               n.cores        <- detectCores()
                               if (opt.estim$Nchains==1)
                               {
                                 self$output  <- self$calibration()
                                 self$mcmc    <- as.mcmc(self$output$out$THETA)
                               } else
                               {
                                 n            <- c(1:opt.estim$Nchains)
                                 self$output  <- mclapply(n,self$calibration,mc.cores = n.cores)
                                 self$mcmc    <- list()
                                 for (i in 1:opt.estim$Nchains)
                                 {
                                   self$mcmc[[i]] <- as.mcmc(self$output[[i]]$out$THETA)
                                 }
                               }
                               cat("\nEnd of the regular calibration\n\n")
                               if (is.null(self$opt.valid)==FALSE)
                               {
                                 cat(paste("\nThe cross validation is currently running on your ",
                                             n.cores," cores available....\n",sep=""))
                                 self$errorCV <- mean(unlist(mclapply(c(1:opt.valid$nCV),
                                                                       self$CV,mc.cores = n.cores)))
                               }
                             },
                             calibration = function(i=NA)
                             {
                               binf          <- private$boundaries()$binf
                               bsup          <- private$boundaries()$bsup
                               MetropolisCpp <- private$MCMC(self$md$model)
                               out           <- MetropolisCpp(self$md$fun,self$opt.estim$Ngibbs,self$opt.estim$Nmh,
                                                              self$opt.estim$thetaInit,self$opt.estim$k,
                                                              self$opt.estim$sig,self$md$Yexp,binf,bsup,self$logPost,1)
                               MAP           <- private$MAPestimator(out$THETA)
                               return(list(out=out,MAP=MAP))
                             },
                             CV = function(i=NA)
                             {
                               if (opt.valid$type.valid=="loo")
                               {
                                 inc <- sample(c(1:self$md$n),1)
                                 if (is.matrix(self$md$X))
                                 {
                                   dataCal <- self$md$X[-inc,]
                                   dataVal <- as.matrix(t(self$md$X[inc,]))
                                 } else
                                 {
                                   dataCal <- self$md$X[-inc]
                                   dataVal <- self$md$X[inc]
                                 }
                                 Ycal        <- self$md$Yexp[-inc]
                                 Yval        <- self$md$Yexp[inc]
                                 mdTemp      <- model(code=self$md$code,dataCal,Ycal,self$md$model,self$md$opt.emul)
                                 mdTempfit   <- self$calibrationCV(mdTemp,Ycal)
                                 Dim         <- length(self$pr)
                                 if (mdTemp$model=="model1" | mdTemp$model=="model2")
                                 {
                                    Ytemp <- mdTemp$pred(mdTempfit$MAP[1:(Dim-1)],mdTempfit$MAP[Dim],dataVal)$y
                                 } else
                                 {
                                   Ytemp <- mdTemp$pred(mdTempfit$MAP[1:(Dim-3)],mdTempfit$MAP[(Dim-3):(Dim-1)]
                                                        ,mdTempfit$MAP[Dim],dataVal)$y
                                 }
                                 return(sqrt((Ytemp-Yval)^2))
                               }
                             },
                             calibrationCV = function(mdTemp,y)
                             {
                               binf          <- private$boundaries()$binf
                               bsup          <- private$boundaries()$bsup
                               MetropolisCpp <- private$MCMC(mdTemp$model)
                               out           <- MetropolisCpp(mdTemp$fun,self$opt.estim$Ngibbs,
                                                              self$opt.estim$Nmh,
                                                              self$opt.estim$thetaInit,self$opt.estim$k,
                                                              self$opt.estim$sig,y,binf,bsup,self$logPost,0)
                               MAP           <- private$MAPestimator(out$THETA)
                               return(list(out=out,MAP=MAP))
                             }
                           ))


calibrate.class$set("private","MCMC",
                  function(model)
                  {
                    switch(model,
                           model1={return(MetropolisHastingsCpp)},
                           model2={return(MetropolisHastingsCpp)},
                           model3={return(MetropolisHastingsCppD)},
                           model4={return(MetropolisHastingsCppD)}
                    )
                  })


calibrate.class$set("private","boundaries",
                function()
                {
                  binf <- self$pr[[1]]$binf
                  bsup <- self$pr[[1]]$bsup
                  for (i in 2:length(self$pr))
                  {
                    binf <- c(binf,self$pr[[i]]$binf)
                    bsup <- c(bsup,self$pr[[i]]$bsup)
                  }
                  return(list(binf=binf,bsup=bsup))
                })


calibrate.class$set("private","logLikelihood",
                function(model)
                {
                  switch(model,
                         model1={return(private$logTest)},
                         model2={return(private$logTest)},
                         model3={return(private$logTestD)},
                         model4={return(private$logTestD)}
                  )
                }
)

calibrate.class$set("private","logTest",
                    function(theta,sig2)
                    {
                      if (length(self$pr) == 1)
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

calibrate.class$set("private","logTestD",
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

calibrate.class$set("private","MAPestimator",
                    function(chain)
                      {
                        dens <- apply(chain,2,density)
                        map <- function(dens)
                        {
                          dens$x[which(dens$y==max(dens$y))]
                        }
                        return(unlist(lapply(dens,map)))
                    })


calibrate.class$set("public","plot",
                    function(graph=c("acf","chains","densities"))
                    {
                      n <- length(self$pr)
                      gg <- list()
                      a <- list()
                      m <- list()
                      p <- list()
                      if ("acf" %in% graph)
                      {
                        for (i in 1:n)
                        {
                          a[[i]] <- self$acf(i)
                        }
                        do.call(grid.arrange,a)
                        gg$acf <- a
                      }
                      if ("chains" %in% graph)
                      {
                        for (i in 1:n)
                        {
                          m[[i]] <- self$mcmcChains(i)
                        }
                        do.call(grid.arrange,m)
                        gg$mcmc <- m
                      }
                      if ("densities" %in% graph)
                      {
                        for (i in 1:n)
                        {
                          p[[i]] <- self$densities(i)
                        }
                        do.call(grid.arrange,p)
                        gg$dens <- p
                      }
                      return(gg)
                    })


calibrate.class$set("public","acf",
                    function(i)
                    {
                      bacf   <- acf(self$output$out$THETA[-c(1:self$opt.estim$burnIn),i], plot = FALSE)
                      bacfdf <- with(bacf, data.frame(lag, acf))
                      p      <- ggplot(data = bacfdf, mapping = aes(x = lag, y = acf))+
                        geom_hline(aes(yintercept = 0))+
                        geom_segment(mapping = aes(xend = lag, yend = 0))+
                        xlab("")+ylab("")+theme_light()
                      return(p)
                      })

calibrate.class$set("public","mcmcChains",
                    function(i)
                    {
                      n <- length(self$output$out$THETA[-c(1:self$opt.estim$burnIn),i])
                      resgg <- data.frame(inc=c(1:n),data=self$output$out$THETA[-c(1:self$opt.estim$burnIn),i])
                      p   <- ggplot(data=resgg, aes(x=inc,y=data))+geom_line()+ylab("")+
                        xlab("")+theme_light()
                      return(p)
                    }
                    )


calibrate.class$set("public","dens",
                    function(i)
                    {
                      binf          <- private$boundaries()$binf
                      bsup          <- private$boundaries()$bsup
                      dplot  <- data.frame(data=self$output$out$THETA[-c(1:self$opt.estim$burnIn),i],type="posterior")
                      p <- ggplot(dplot,aes(data,fill=type,color=type)) +
                        geom_density(kernel = "gaussian",adjust=3,alpha=0.1)+
                        theme_light()+xlab("")+ylab("")+ xlim(binf[i],bsup[i])+
                        scale_fill_manual( values = "blue")+
                        scale_color_manual(values = "blue")+
                        theme(legend.position=c(0.86,0.86),
                              legend.text=element_text(face="bold",size = '20'),
                              legend.title=element_blank(),
                              legend.key=element_rect(colour=NA),
                              axis.text=element_text(size=20))+
                        geom_hline(aes(yintercept = 0))
                      return(p)
                    }
)


