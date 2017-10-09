#' A Reference Class to generate different calibration methods after generating a model from
#' model.class
#'
#' @examples
#'
#' @export
estim.class <- R6::R6Class(classname = "estim.class",
                 public = list(
                    code       = NULL,
                    X          = NULL,
                    Yexp       = NULL,
                    model      = NULL,
                    type.prior = NULL,
                    log        = NULL,
                    opt.emul   = NULL,
                    opt.prior  = NULL,
                    opt.estim  = NULL,
                    binf       = NULL,
                    bsup       = NULL,
                    md         = NULL,
                    pr         = NULL,
                    initialize = function(code=NA,X=NA,Yexp=NA,model=NA,type.prior=NA,log=TRUE,
                                          opt.emul=NA,opt.prior=NA,opt.estim=NA)
                    {
                      self$code       <- code
                      self$X          <- X
                      self$Yexp       <- Yexp
                      self$model      <- model
                      self$log        <- log
                      self$opt.emul   <- opt.emul
                      self$opt.prior  <- opt.prior
                      self$opt.estim  <- opt.estim
                      self$type.prior <- type.prior
                      private$boundaries()
                      private$checkup()
                      self$md         <- model(code,X,Yexp,model,opt.emul)
                      self$pr         <- prior(self$type.prior,opt.prior,log=TRUE)
                    },
                    estimation = function()
                    {
                      print(self$logTest(10,1))
                      out <- MetropolisHastingsCpp(self$md$fun,self$opt.estim$Ngibbs,
                                                   self$opt.estim$Nmh,self$opt.estim$thetaInit,
                                                   self$opt.estim$k,self$opt.estim$sig,self$Yexp,
                                                   self$binf,self$bsup,private$logTest)
                      return(out)
                    }
                   ))


estim.class$set("private","boundaries",
                function()
                {
                  self$binf <- self$pr[[1]]$binf
                  self$bsup <- self$pr[[1]]$bsup
                  for (i in 2:length(self$type.prior))
                  {
                    self$binf <- cbind(self$binf,self$pr[[i]]$binf)
                    self$bsup <- cbind(self$bsup,self$pr[[i]]$bsup)
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



















