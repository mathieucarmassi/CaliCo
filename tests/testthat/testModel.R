rm(list=ls())
X <- seq(0,1,length.out=5)
code <- function(X,theta)
{
  return((6*X-2)^2*sin(theta*X-6))
}
Yexp <- code(X,10.5)+rnorm(5,0,0.1)
model1 <- model(code,X,Yexp,"model1")
expect_is(model1,"model.class")
binf <- 10
bsup <- 11
model2 <- model(code,X,Yexp,"model2",
                opt.gp=list(type="matern5_2", DOE=NULL),
                opt.emul=list(p=1,n.emul=10,binf=binf,bsup=bsup,type="maximinLHS"))
expect_is(model2,"model.class")
DOE <- DiceDesign::lhsDesign(20,2)$design
DOE[,2] <- unscale(DOE[,2],binf,bsup)
model2 <- model(code,X,Yexp,"model2",
                opt.gp=list(type="matern5_2", DOE=DOE))
expect_is(model2,"model.class")
Ysim <- matrix(nr=20,nc=1)
for (i in 1:20)
{
  covariates <- as.matrix(DOE[i,1])
  Ysim[i] <- code(covariates,DOE[i,2])
}
model2 <- model(code=NULL,X,Yexp,"model2",
                opt.gp = list(type="matern5_2", DOE=NULL),
                opt.sim=list(Ysim=Ysim,DOEsim=DOE))
expect_is(model2,"model.class")
model3 <- model(code,X,Yexp,"model3",
                opt.disc=list(kernel.type="gauss"))
expect_is(model3,"model.class")
model4 <- model(code=NULL,X,Yexp,"model4",
                opt.disc=list(kernel.type="gauss"),
                opt.gp = list(type="matern5_2", DOE=NULL),
                opt.sim=list(Ysim=Ysim,DOEsim=DOE))
expect_is(model4,"model.class")

