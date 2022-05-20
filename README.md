# MTPLEN
The survival data of outliers were analyzed by robust Cox regression model based on truncation, and the model was fitted by truncation of individuals with small contributions to partial likelihood, so as to obtain robust estimates unaffected by outliers.
# Intallation
The required R package needs to be installed before analysis,You can copy the following code. 
install.packages(c("glmne","survival","mvtnorm","glmnetUtils","glmpath",”MASS”,"pacman")) 
library(pacman) 
p_load(glmnet,survival,mvtnorm,glmnetUtils,glmpath,MASS)
install_github("hwsun2000/MTPLEN")
# usage
fitter(t,delta,x,gamma,poscaso=NULL,D=NULL,kmax=NULL,iter,method,rseqmax=NULL,n.multistart=NULL)
# Example
library(glmnet)
library(survival)
library(mvtnorm)
library(glmnetUtils)
library(glmpath)
library(MASS)

n=100 
p=200
q <- 12
p <- 200
CL<- 1
CU<- 3
epsilon=0.1
beta1=rep(0.5,10)
beta0=rep(0,p-length(beta1))
beta=c(beta1,beta0)

seed<- 8601
set.seed(seed)
ssigma =0.7 Sigma =ssigma^t(sapply(1:p, function(i, j) abs(i-j), 1:p)) x=rmvnorm(n, sigma=Sigma) ht=exp(x%*%beta)
maxht=max(ht)
minht=min(ht)
kk=sample(n,floor(n*epsilon))

seed<- 8007
set.seed(seed)
flag=rbinom(length(kk),1,0.5)
ht[kk]=minht*flag+(maxht)*(1-flag)
routliers=kk
rinliers=(1:n)[-routliers]

#observed survivial time
logS=log(matrix(runif(n,0,1),n,1)) #log[S(t)]
T=-logS/ht  #survival time

#censored time
myrate <- matrix(runif(n,CL,CU),n,1)/ht
C <- rexp(n, rate = 1/myrate)

#survival time and state
time <- apply(cbind(C,T),1,min)
t<-time
stat <- (T<C)
surv=Surv(time,stat)

mtlcox<-fitter(time,stat,x,gamma=0.25,iter="EN",method="EN",n.multistart=1)
mtlcox$rwt.betas  #Coefficients of RwtMTPLEN
mtlcox$rwt.outliers  #Outliers or influential observations
mtlcox$devres    #deviance residual of observaiton
mtlcox$rwt.max.lik   #likelihood function of RwtMTPLEN
# Development
This R package are develpoed by Hongwei Sun and wenjing Zhang.
