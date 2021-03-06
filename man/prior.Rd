% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{prior}
\alias{prior}
\title{Generates \code{\link{prior.class}} objects.}
\usage{
prior(type.prior, opt.prior)
}
\arguments{
\item{type.prior}{the vector of the prior types selected. For example type.prior=c("gaussian","gamma")}

\item{opt.prior}{list of the hyperparameters relatives to the prior selected. If the first prior selected is
Gaussian, the hyperparameters would be the mean and the standard deviation. See Details for further clarifications.}
}
\value{
\code{prior} returns a \code{\link{prior.class}} object. Two main methods are available:
\itemize{\item{plot()}{ display the probability density of the prior}
\item{print()}{ returns the main information concerning the prior distribution}}
}
\description{
\code{prior} is a function that generates a \code{\link{prior.class}} containing information about one or
several priors. When several priors are selected, the function \code{prior}
 returns a list of \code{\link{prior.class}}.
}
\details{
The densities implemented are defined as follow
\itemize{
\item{The Gaussian density:
\deqn{f(x)=1/(\sigma*\sqrt(2\pi))exp{-1/2*((x-\mu)/\sigma)^2}}
where \strong{\eqn{\mu}} and \strong{\eqn{\sigma}} (the mean and the standard deviation)
are the two hyperparameters. The vector \eqn{c(\mu,\sigma^2)} is the one looked for in opt.prior.}
\item{The Gamma density:
\deqn{f(x)=1/(k^a*\Gamma(a))*x^(a-1)*exp(-(x/k))}
where \strong{\eqn{a}} and \strong{\eqn{k}} (the shape and the scale)
are the two hyperparameters. The vector \eqn{c(a,k)} is the one looked for in opt.prior.}
\item{The Uniform density:
\deqn{f(x)=1/(b-a)}
where \strong{\eqn{a}} and \strong{\eqn{b}} (the upper and the lower bound)
are the two hyperparameters. The vector \eqn{c(a,b)} is the one looked for in opt.prior.}
}
}
\examples{
\dontrun{
#### Only one prior is wanted
## For a Gaussian Prior
gaussian <- prior(type.prior="gaussian",opt.prior=list(c(0.5,0.001)))
plot(gaussian)

#### For several priors
priors <- prior(type.prior=c("gaussian","gamma"),opt.prior=list(c(0.5,0.001),c(5,1)))
plot(priors$Prior1)
plot(priors$Prior2)
}
}
\seealso{
\code{\link{model}}, \code{\link{calibrate}}, \code{\link{forecast}}, \code{\link{sequentialDesign}}
}
\author{
M. Carmassi
}
