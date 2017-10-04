#' A Reference Class to generates differents prior objects
#'
#' @description See the function blabla which produces an instance of this class
#' This class comes with a set of methods, some of them being useful for the user:
#' See the documentation for blabla... Other methods
#'  should not be called as they are designed to be used during the optimization process.
#'
#' Fields should not be changed or manipulated by the user as they are updated internally
#' during the estimation process.
#' @field type.prior of the selected prior
#' @field opt.prior the charasteristics of the selected prior
#' @field log if we want the log result or not
prior.class <- R6::R6Class(classname = "prior.class",
                           public = list(
                             type.prior  = NULL,
                             opt.prior   = NULL,
                             log         = NULL,
                             initialize = function(type.prior=NA,opt.prior=NA,log=NA)
                             {
                               self$type.prior  <- type.prior
                               self$opt.prior   <- opt.prior
                               self$log         <- log
                               private$checkPrior()
                               private$loadPackages()
                             }
                           ))

prior.class$set("private","checkPrior",
                function()
                {
                  if (self$type.prior != "gaussian" && self$type.prior != "gamma" &&
                      self$type.prior != "invGamma" && self$type.prior != "unif")
                  {
                    stop('Please select a correct prior')
                  }
                })

prior.class$set("private","loadPackages",
                function()
                {
                  library(R6)
                  library(DiceDesign)
                  library(DiceKriging)
                  library(FactoMineR)
                })


gaussian.class <- R6::R6Class(classname = "gaussian.class",
                              inherit = prior.class,
                              public = list(
                                mean  = NULL,
                                var   = NULL,
                                initialize = function(type.prior=NA,opt.prior=NA,log=NA)
                                {
                                  super$initialize(type.prior,opt.prior,log)
                                  self$mean  <- unlist(self$opt.prior)[1]
                                  self$var   <- unlist(self$opt.prior)[2]
                                },
                                gaussian = function(x=seq(0,1,length.out = 100))
                                  {
                                  if (self$log == FALSE)
                                  {
                                    return(1/(sqrt(2*pi*self$var))*
                                               exp(-1/(2*self$var)*
                                                     (x-self$mean)^2))
                                    } else
                                  {
                                    return(log(1/(sqrt(2*pi*self$var))*
                                                   exp(-1/(2*self$var)*
                                                         (x-self$mean)^2)))
                                  }
                                }
                              ))


unif.class <- R6::R6Class(classname = "unif.class",
                              inherit = prior.class,
                              public = list(
                                binf    = NULL,
                                bsup    = NULL,
                                y       = NULL,
                                initialize = function(type.prior=NA,opt.prior=NA,log=NA)
                                {
                                  super$initialize(type.prior,opt.prior,log)
                                  self$binf <- unlist(self$opt.prior)[1]
                                  self$bsup <- unlist(self$opt.prior)[2]
                                },
                                unif = function(x=seq(0,1,length.out = 100))
                                {
                                  self$y <- matrix(nr=length(x),1)
                                  for (i in 1:length(x))
                                  {
                                    if (x[i]>=self$binf & x[i]<=self$bsup)
                                    {
                                      self$y[i] <- 1/(self$bsup-self$binf)
                                    } else
                                    {
                                      self$y[i] <- 0
                                    }
                                  }
                                  if (self$log == FALSE)
                                  {
                                    return(self$y)
                                  } else
                                  {
                                    return(log(self$y))
                                  }
                                }
                              ))


gamma.class <- R6::R6Class(classname = "gamma.class",
                          inherit = prior.class,
                          public = list(
                            shape    = NULL,
                            scale    = NULL,
                            y       = NULL,
                            initialize = function(type.prior=NA,opt.prior=NA,log=NA)
                            {
                              super$initialize(type.prior,opt.prior,log)
                              self$shape <- unlist(self$opt.prior)[1]
                              self$scale <- unlist(self$opt.prior)[2]
                            },
                            gamma = function(x=seq(0,1,length.out = 100))
                            {
                              if (self$log == FALSE)
                              {
                                return(1/(self$scale^self$shape*gamma(self$shape))*
                                         x^(self$shape-1)*exp(-(x/self$scale)))
                              } else
                              {
                                return(log(1/(self$scale^self$shape*gamma(self$shape))*
                                             x^(self$shape-1)*exp(-(x/self$scale))))
                              }
                            }
                          ))


invGamma.class <- R6::R6Class(classname = "invGamma.class",
                           inherit = prior.class,
                           public = list(
                             shape    = NULL,
                             scale    = NULL,
                             y       = NULL,
                             initialize = function(type.prior=NA,opt.prior=NA,log=NA)
                             {
                               super$initialize(type.prior,opt.prior,log)
                               self$shape <- unlist(opt.prior)[1]
                               self$scale <- unlist(opt.prior)[2]
                             },
                             invGamma = function(x=seq(0,1,length.out = 100))
                             {
                               if (self$log == FALSE)
                               {
                                 return(1/(1/(self$scale^self$shape*gamma(self$shape))*
                                          x^(self$shape-1)*exp(-(x/self$scale))))
                               } else
                               {
                                 return(log(1/(1/(self$scale^self$shape*gamma(self$shape))*
                                              x^(self$shape-1)*exp(-(x/self$scale)))))
                               }
                             }
                           ))


