## ---- fig.width=4,echo=FALSE,fig.align='center'--------------------------
library(ggplot2)
n=20
X <- seq(0,1,length.out=n)
code <- function(X,theta)
{
  return((6*X-2)^2*sin(theta*X-4))
}
Yexp <- code(X,10.5)+rnorm(n,0,sqrt(0.01))
Yt <- code(X,11)
gdata <- data.frame(y=Yexp,x=X,type="Experiments")
gdata2 <- data.frame(y=Yt,x=X,type="Theoretical results")
gdata <- rbind(gdata,gdata2)

p1 <- ggplot(gdata,aes(x = x,y=y, col=type))+geom_line() + theme_light() + theme(legend.position=c(0.50,0.80),
                                                                                 legend.title=element_blank())
p1

## ---- fig.width=6,echo=FALSE,fig.align='center'--------------------------
theta.prior <- rnorm(100,11,sqrt(.5))
sigma.prior <- 0.01

Y_prior <- matrix(nr=100,nc=n)
for (i in 1:100){
  Y_prior[i,] <- code(X,theta.prior[i])+rnorm(n,0,sqrt(sigma.prior))
}

qq <- apply(Y_prior,2,quantile,probs=c(0.05,0.95))
gdata <- data.frame(y=Yexp,x=X,upper=qq[2,],lower=qq[1,],type="Experiments",fill="credibility interval at 90%")
gdata2 <- data.frame(y=Yt,x=X,upper=qq[2,],lower=qq[1,],type="Prior belief",fill="credibility interval at 90%")
gdata <- rbind(gdata,gdata2)

p2 <- ggplot(gdata)+geom_line(aes(x=x,y=y,col=type))+geom_ribbon(aes(x=x,ymax=upper,ymin=lower,fill=fill),alpha=0.3)+theme_light()+theme(legend.title=element_blank(),legend.key=element_rect(colour=NA),legend.text=element_text(size = '8'))
p2

## ---- echo=TRUE----------------------------------------------------------
library(CaliCo)
# Number of experiments
n=20
# Time interval
t <- seq(0,1,length.out=n)
# Code definition
code <- function(t,theta)
{
  return((6*t-2)^2*sin(theta*t-4))
}
# Generate the experiment
Yexp <- code(t,10.5)+rnorm(n,0,sqrt(0.01))
# Generate the first model
model1 <- model(code=code,X=t,Yexp=Yexp,model="model1")

## ----echo=TRUE-----------------------------------------------------------
model1$print()

## ---- echo=TRUE,fig.width=4,fig.align="center"---------------------------
model1$plot(11,0.01)

## ---- echo=TRUE,fig.width=4,fig.align="center"---------------------------
p <- model1$plot(11,0.01)
p+ggtitle("Model1 and experiments")+ylab("y")+xlab("t")+
  theme(legend.position=c(0.6,0.75),legend.text=element_text(size = '7'))

## ----echo=TRUE-----------------------------------------------------------
# Set the lower and upper bound (taken large intentionally)
binf <- 8
bsup <- 13

# Set the emulation option
opt.emul <- list(p=1,n.emul=30,type="matern3_2",binf=binf,bsup=bsup,DOE=NULL)
# Generate the second model
model2 <- model(code,X,Yexp,"model2",opt.emul)

## ---- echo=TRUE,fig.width=4,fig.align="center"---------------------------
p <- model2$plot(11,0.01,points=TRUE)
p+ylab("y")+xlab("t")+
  theme(legend.position=c(0.6,0.75),legend.text=element_text(size = '7'))

## ---- echo=TRUE----------------------------------------------------------
opt.disc <- list(kernel.type="matern5_2")
# Generate the third model
model3 <- model(code,X,Yexp,"model3",opt.disc=opt.disc)

# Generate the fourth model where the emulation option are needed
### Here opt.emul is the same as the one for the second model
model4 <- model(code,X,Yexp,"model4",opt.emul,opt.disc)

## ---- echo=TRUE, fig.width=4, fig.align="center"-------------------------
p <- model3$plot(11,c(2,0.5),0.01)
p+ylab("y")+xlab("t")+
  theme(legend.position=c(0.6,0.75),legend.text=element_text(size = '7'))

## ---- echo=TRUE, fig.width=6, fig.align="center"-------------------------
p <- model4$plot(11,c(1,0.1),0.01)
p+ylab("y")+xlab("t")+
  theme(legend.position="right",legend.text=element_text(size = '7'))

## ---- fig.width=4,fig.align='center'-------------------------------------
gaussian <- prior(type.prior="gaussian",opt.prior=list(c(0.5,0.001)))
p <- gaussian$plot()
p

## ---- fig.width=4,fig.height=5,fig.align='center'------------------------
priors <- prior(type.prior=c("gaussian","gamma"),opt.prior=list(c(0.5,0.001),c(5,1)))
grid.arrange(priors$Prior1$plot(),priors$Prior2$plot(),nrow=2)

## ---- echo=TRUE----------------------------------------------------------
priors$Prior1
priors$Prior2

