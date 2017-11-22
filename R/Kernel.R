#' A Reference Class to generates differents model objects
#'
#' @description See the function \code{\link{model}} which produces an instance of this class
#'
#' @export
Kernel.class <- R6::R6Class(classname = "Kernel.class",
                           public = list(
                             X           = NULL,
                             var         = NULL,
                             theta       = NULL,
                             n           = NULL,
                             Kernel.type = NULL,
                             Cov         = NULL,
                             initialize = function(X=NA,var=NA,theta=NA,Kernel.type=NA)
                             {
                               self$X           <- X
                               self$var         <- var
                               self$theta       <- theta
                               if (is.null(dim(self$X))==TRUE)
                               {
                                 self$n <- length(self$X)
                               } else
                               {
                                 self$n <- nrow(self$X)
                               }
                               self$Kernel.type <- Kernel.type
                             }
                           ))


gauss.class <- R6::R6Class(classname= "gauss.class",
                           inherit = Kernel.class,
                           public = list(
                             initialize = function(X,var,theta,Kernel.type="gauss")
                             {
                               super$initialize(X,var,theta,Kernel.type="gauss")
                               self$Cov <- self$covariance()
                             },
                             covariance = function()
                             {
                               self$Cov <- matrix(nr=self$n,nc=self$n)
                               if (is.null(ncol(self$X)))
                               {
                                 self$Cov <- matrix(nr=self$n,nc=self$n)
                                 for (j in 1:self$n)
                                 {
                                   for (i in 1:self$n)
                                   {
                                     self$Cov[i,j] <- self$var*exp(-1/2*(sum((self$X[i]-self$X[j])^2)/self$theta)^2)
                                   }
                                 }
                               } else
                               {
                                 self$Cov <- matrix(nr=self$n,nc=self$n)
                                 for (j in 1:self$n)
                                 {
                                   for (i in 1:self$n)
                                   {
                                     self$Cov[i,j] <- self$var*exp(-1/2*(sum((self$X[i,]-self$X[j,])^2)/self$theta)^2)
                                   }
                                 }
                               }
                               return(self$Cov)
                             }
                           ))

exp.class <- R6::R6Class(classname= "exp.class",
                           inherit = Kernel.class,
                         public = list(
                           initialize = function(X,var,theta,Kernel.type="exp")
                           {
                             super$initialize(X,var,theta,Kernel.type="exp")
                             self$Cov <- self$covariance()
                           },
                           covariance = function(X,var,theta)
                           {
                             self$Cov <- matrix(nr=self$n,nc=self$n)
                             if (is.null(ncol(self$X)))
                             {
                               self$Cov <- matrix(nr=self$n,nc=self$n)
                               for (j in 1:self$n)
                               {
                                 for (i in 1:self$n)
                                 {
                                   self$Cov[i,j] <- self$var*exp(sum((self$X[i]-self$X[j])^2)/self$theta)
                                 }
                               }
                             } else
                             {
                               self$Cov <- matrix(nr=self$n,nc=self$n)
                               for (j in 1:self$n)
                               {
                                 for (i in 1:self$n)
                                 {
                                   self$Cov[i,j] <- self$var*exp(sum((self$X[i,]-self$X[j,])^2)/self$theta)
                                 }
                               }
                             }
                             return(self$Cov)
                           }
))


matern3_2.class <- R6::R6Class(classname= "matern3_2.class",
                         inherit = Kernel.class,
                         public = list(
                         initialize = function(X,var,theta,Kernel.type="matern3_2")
                         {
                           super$initialize(X,var,theta,Kernel.type="matern3_2")
                           self$Cov <- self$covariance()
                         },
                         covariance = function(X,var,theta)
                         {
                           self$Cov <- matrix(nr=self$n,nc=self$n)
                           if (is.null(ncol(self$X)))
                           {
                             self$Cov <- matrix(nr=self$n,nc=self$n)
                             for (j in 1:self$n)
                             {
                               for (i in 1:self$n)
                               {
                                 self$Cov[i,j] <- self$var*(1+sqrt(3)*sum((self$X[i]-self$X[j])^2)/self$theta)*
                                   exp(-sqrt(3)*sum((self$X[i]-self$X[j])^2)/self$theta)
                               }
                             }
                           } else
                           {
                             self$Cov <- matrix(nr=self$n,nc=self$n)
                             for (j in 1:self$n)
                             {
                               for (i in 1:self$n)
                               {
                                 self$Cov[i,j] <- self$var*(1+sqrt(3)*sum((self$X[i,]-self$X[j,])^2)/self$theta)*
                                   exp(-sqrt(3)*sum((self$X[i,]-self$X[j,])^2)/self$theta)
                               }
                             }
                           }
                           return(self$Cov)
                         }
))


matern5_2.class <- R6::R6Class(classname= "matern5_2.class",
                               inherit = Kernel.class,
                               public = list(
                               initialize = function(X,var,theta,Kernel.type="matern5_2")
                               {
                                 super$initialize(X,var,theta,Kernel.type="matern5_2")
                                 self$Cov <- self$covariance()
                               },
                               covariance = function(X,var,theta)
                               {
                                 self$Cov <- matrix(nr=self$n,nc=self$n)
                                 if (is.null(ncol(self$X)))
                                 {
                                   self$Cov <- matrix(nr=self$n,nc=self$n)
                                   for (j in 1:self$n)
                                   {
                                     for (i in 1:self$n)
                                     {
                                       self$Cov[i,j] <- self$var*(1+sqrt(5)*sum((self$X[i]-self$X[j])^2)/self$theta+
                                                               5/3*(sum((self$X[i]-self$X[j])^2)/self$theta)^2)*
                                                               exp(-sqrt(5)*sum((self$X[i]-self$X[j])^2)/self$theta)
                                     }
                                   }
                                 } else
                                 {
                                   self$Cov <- matrix(nr=self$n,nc=self$n)
                                   for (j in 1:self$n)
                                   {
                                     for (i in 1:self$n)
                                     {
                                       self$Cov[i,j] <- self$var*(1+sqrt(5)*sum((self$X[i,]-self$X[j,])^2)/self$theta+
                                                               5/3*(sum((self$X[i,]-self$X[j,])^2)/self$theta)^2)*
                                                               exp(-sqrt(5)*sum((self$X[i,]-self$X[j,])^2)/self$theta)
                                     }
                                   }
                                 }
                                 return(self$Cov)
                               }
))
