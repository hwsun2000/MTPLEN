
#' Robust penalized Cox regression based on truncation
#'
#' @param t the time variables.
#' @param delta event indicators.
#' @param x matrix of covariates.
#' @param gamma trimming level.
#' @param poscaso optional initial solution.
#' @param D tuning parameters for the fitting algorithm.
#' @param kmax tuning parameters for the fitting algorithm.
#' @param iter the penalized method used in iteration.There are two selection:"EN" or "LASSO".
#' @param method the penalized method used for last optimal subset.There are two selection:"EN" or "LASSO".
#' @param rseqmax tuning parameters for the fitting algorithm.
#' @param n.multistart tuning parameters for the fitting algorithm.
#'
#' @return max.beta, link_raw, rwt.betas, devres, rwt.indices, rwt.max.lik, link_rwtb, link_raw

library(survival)
library(glmpath)
library(SIS)
library(glmnetUtils)
library(MASS)

fitter=function(t,delta,x,gamma,poscaso=NULL,D=NULL,kmax=NULL,iter,method,rseqmax=NULL,n.multistart=NULL) {

# t: times
# delta: event indicators
# x: matrix of covariates
# poscaso: optional initial solution
# gamma: trimming level
# D, kmax, rseqmax, n.multistart: tuning parameters for the fitting algorithm
# iter, the penalized method used in iteration.There are two selection:"EN" or "LASSO"
# method, the penalized method used for last optimal subset.There are two selection:"EN" or "LASSO".

n=length(t)
p=dim(x)[2]


hn=floor((n+1)*(1-gamma))

if(is.null(kmax)) {kmax=1000}
if(!is.matrix(x)) {x=matrix(x)}
if(is.null(rseqmax)) {rseqmax=1}
if(is.null(D)) {D=0.1}
if(is.null(poscaso)) {poscaso=sample(n,ceiling(n*(1-gamma)))}
if(is.null(n.multistart)) {n.multistart=1}
x=as.matrix(x)
formersubset=500
ncsteps=2
latersubset=10


max.lik=-10000

n1=length(which(delta==1))
n0=length(which(delta==0))
index1<-which(delta==1)
index0<-which(delta==0)

sn1=floor((n1+1)*(1-gamma))
sn0=hn-sn1

#normalize by  mean and standard deviation

       center=apply(x,2,mean)
       x.c = sweep(x, 2, center)
       sd.c=apply(x,2,sd)
       xsd=sweep(x.c,2,sd.c,"/")


   margcoef = abs(cor(x,t))
   maxcc=max(margcoef)
   frac=seq(0.05,1,by=0.05)*maxcc  # 20??lambda


##????logplik

logplik <- function(x, time, status, b, method = c('breslow', 'efron'),
                    return.all = FALSE)
  {
    method <- match.arg(method)
    n <- length(time)
    o <- order(status, decreasing=T)
    oo <- o[order(time[o])]
    time <- time[oo]
    status <- status[oo]
    rept <- rep(0, n)
    for (i in 1:n) rept[i] <- sum(time[i:n] == time[i] & status[i:n] == 1)
    complete <- which(status == 1)
    nnc <- length(complete)
    if (nnc == 0) {
      stop('No complete observation. Failed to compute partial likelihood.')
    }
    dmat <- matrix(0, n, nnc)
    for (i in 1:nnc) {
      dmat[time >= time[complete[i]], i] <- 1
      if (method == 'efron') {
        if (rept[complete[i]] > 0) {
          tie <- time == time[complete[i]] & status == 1
          di <- max(rept[tie])
          dmat[tie, i] <- dmat[tie, i] - (di - rept[complete[i]]) / di
        }
      }
    }
    eta <- x %*% b
    eeta <- exp(eta)
    k <- ncol(eta)
    loglik <- rep(0, k)
    for (i in 1:k) {
      w <- dmat * eeta[oo, i]
      wsum <- apply(w, 2, sum)
      loglik[i] <- sum(eta[oo, i][status == 1]) - sum(log(wsum))
    }
    if (return.all) {
      return(list(loglik = loglik, w = scale(w, F, wsum), eta = eta,
                  dmat = dmat, oo = oo))
    } else {
      return(loglik)
    }
  }


   ###2 Csteps######

      cstep<-function(indices){

     for (i in 1:ncsteps){

     #normalize by trimmed mean and standard

      center.trim=apply(x[indices,],2,mean)    ##trimmed mean and SD
      sd.trim=apply(x[indices,],2,sd)

    sd.tc<-ifelse(sd.trim==0,sd.c,sd.trim)
    center.tc<-ifelse(sd.trim==0,center,center.trim)

       x.c = sweep(x, 2, center.tc)
       xsd=sweep(x.c,2,sd.tc,"/")

        xsd.trim=xsd[indices,]
        y.trim=Surv(t[indices],delta[indices])

         ###################Deviance residual###############
        devres=rep(1,n)

        stime=cv.glmnet(x=xsd.trim,y=y.trim,family='cox',alpha =1,standardize=FALSE)
        trim.beta=as.numeric(coef(stime,s="lambda.1se"))
        link_EN=as.matrix(predict(stime,newx=xsd.trim,type="link"))

        strata=rep(1,length(indices))
        time <-t[indices]
        status <- delta[indices]
        sorted <- order(strata, time)
        strata <- strata[sorted]
	newstrat <- as.integer(c(1*(diff(as.numeric(strata))!=0), 1))
        stime <- as.double(time[sorted])
        sstat <- as.integer(status[sorted])
        ### method=efron ##
        score=exp(link_EN)[sorted]
        weights=rep(1,length(indices))
        coxres <- .C("coxmart", as.integer(length(indices)),
				1,
				stime,
				sstat,
				newstrat,
				as.double(score),
				as.double(weights),
				resid=double(length(indices)))
            resid <- double(length(indices))
            resid[sorted] <- coxres$resid     ##Matingaile residual

          	devres[indices]=sign(resid) *sqrt(-2* (resid+
			     ifelse(status==0, 0, status*log(status-resid))))

           score_all=exp(xsd%*%trim.beta)
           score[score<=0.0000001]=0.0000001
           h0=(status-resid)/exp(link_EN)    ###baseline risk
            h0[h0<=0.0000001]=0.0000001

           hh<-as.data.frame(cbind(time,h0))

           myfit<-lqs(log(h0)~log(time),data=hh,method="lms")

           h00=exp(as.vector(coefficients(myfit))[1]+as.vector(coefficients(myfit))[2]*log(t[-indices]))

           mres=delta[-indices]-h00*score_all[-indices]

          devres[-indices]=sign(mres) *sqrt(-2* (mres+
			     ifelse(delta[-indices]==0, 0, delta[-indices]*log(delta[-indices]-mres))))

           o1<-index1[order(abs(devres[index1]),decreasing=FALSE)[1:sn1]]

           o0<-index0[order(abs(devres[index0]),decreasing=FALSE)[1:sn0]]

         indices<-sort(c(o1,o0))

     }

     ### for subset, applied penalized COX regression

         #normalize by trimmed mean and standard

      center.trim=apply(x[indices,],2,mean)
      sd.trim=apply(x[indices,],2,sd)

      sd.tc<-ifelse(sd.trim==0,sd.c,sd.trim)
      center.tc<-ifelse(sd.trim==0,center,center.trim)

       x.c = sweep(x, 2, center.tc)
       xsd=sweep(x.c,2,sd.tc,"/")


       xsd.trim=xsd[indices,]
       y.trim=Surv(t[indices],delta[indices])

       stime=cv.glmnet(x=xsd.trim,y=y.trim,family='cox',alpha =1,standardize=FALSE)

       trim.beta=as.matrix(coef(stime,s="lambda.1se"))

       link_EN=as.matrix(predict(stime,newx=xsd.trim,type="link"))

       crit=logplik(xsd.trim,t[indices],delta[indices],trim.beta)

     result<-list(indices=indices,crit=crit,mybetas=trim.beta)
     return(result)

  }


  #### Draw subsets and perform 2 AR-Csteps#####
   subset2=matrix(0,20,formersubset)

   for (j in 1:formersubset)
  {
   subset2[,j]=sample(1:n,20)
   }

   statn<-as.numeric(stat)

   statmat<-matrix(0,20,formersubset)

     for (j in 1:formersubset)
  {
   statmat[,j]=statn[subset2[,j]]
   }
    statsum<- apply(statmat, 2, sum)
    subset2<- subset2[,which(statsum>2)]

   multi_indices=matrix(0,hn,ncol(subset2))
   multi_crit=rep(0,ncol(subset2))
   xishu=matrix(0,p,ncol(subset2))

   # perform nsctep C-step on formersubset samples
   for (k in 1:ncol(subset2))
   {
    jieguo<-cstep(subset2[,k])
    multi_crit[k]<-jieguo$crit    #likelihood function
    multi_indices[,k]<-jieguo$indices  #subscript of subsamples
    xishu[,k]<-jieguo$mybetas  # coefficient
    }



  #keep latersubset the optimal subset with largest likelihood

     betterindices_indi=order(multi_crit,decreasing=TRUE)[1:latersubset]
     betterindices=multi_indices[,betterindices_indi]


  cstep2<-function(indices){

     poscaso=indices

    t.sub=t[poscaso]
    delta.sub=delta[poscaso]
       y.sub = survival::Surv(t.sub,delta.sub)

      #normalize by trimmed mean and standard

      center.trim=apply(x[poscaso,],2,mean)
      sd.trim=apply(x[poscaso,],2,sd)

    sd.tc<-ifelse(sd.trim==0,sd.c,sd.trim)
    center.tc<-ifelse(sd.trim==0,center,center.trim)

       x.c = sweep(x, 2, center.tc)
       xsd=sweep(x.c,2,sd.tc,"/")

        xsd.sub=xsd[poscaso,]


     if(iter == "EN")
    {
     my.alpha <- seq(0.1,1,0.1)
     cvm=rep(0,length(my.alpha))
     var.selected.EN <- matrix(0,dim(xsd)[2],length(my.alpha))
     link<- matrix(0,length(poscaso),length(my.alpha))

     fit.EN.cv=cva.glmnet(x=xsd.sub,y=y.sub,alpha=my.alpha,lambda=frac,family="cox",standardize=FALSE)

      for (j in 1:length(my.alpha)){

           cvm[j]=fit.EN.cv$modlist[[j]]$cvm[which(fit.EN.cv$modlist[[j]]$cvm==min(fit.EN.cv$modlist[[j]]$cvm))]

          var.selected.EN[,j]=as.matrix(coef(fit.EN.cv$modlist[[j]],s="lambda.min"))
       link[,j]=xsd.sub%*%var.selected.EN[,j]
       }


      beta=var.selected.EN[,which(cvm==min(cvm))]
       link_EN=link[,which(cvm==min(cvm))]


      } else if (iter == "LASSO") {

      stime=cv.glmnet(x=xsd.sub,y=y.sub,family='cox',alpha =0.3,standardize=FALSE)

       beta=coef(stime,s="lambda.1se")
     link_EN=as.matrix(predict(stime,newx=xsd.sub,type="link"))

        }



#??DDlogpliko????y
lik=logplik(xsd.sub,t.sub,delta.sub,beta)

likold=-10000
k=1
bestlik=lik
b2=beta
p2=poscaso
rseq=0

while(k<kmax && rseq<rseqmax) {


Tk=log(k+1,2)/D
k=k+1




if(lik-likold!=0) {rseq=0}
if(lik-likold==0) {rseq=rseq+1}

likold=lik

  ####Baseline risks are estimated by subsets and the residuals of all individuals are estimated accordingly##################

     devres=rep(0,n)

     strata=rep(1,length(poscaso))
     time <-t[poscaso]
     status <- delta[poscaso]
     sorted <- order(strata, time)
	strata <- strata[sorted]
	newstrat <- as.integer(c(1*(diff(as.numeric(strata))!=0), 1))
      stime <- as.double(time[sorted])
      sstat <- as.integer(status[sorted])

      score=exp(link_EN)[sorted]
      weights=rep(1,length(poscaso))
      coxres <- .C("coxmart", as.integer(length(poscaso)),
				1,
				stime,
				sstat,
				newstrat,
				as.double(score),
				as.double(weights),
				resid=double(length(poscaso)))
            resid <- double(length(poscaso))
            resid[sorted] <- coxres$resid     ##Matingaile residual

          	devres[poscaso]=sign(resid) *sqrt(-2* (resid+
			     ifelse(status==0, 0, status*log(status-resid))))

           score_all=exp(xsd%*%beta)
           score[score<=0.0000001]=0.0000001
           h0=(status-resid)/exp(link_EN)
            h0[h0<=0.0000001]=0.0000001

            hh<-as.data.frame(cbind(time,h0))
           myfit<-lqs(log(h0)~log(time),data=hh,method="lms")
           h00=exp(as.vector(coefficients(myfit))[1]+as.vector(coefficients(myfit))[2]*log(t[-poscaso]))
           mres=delta[-poscaso]-h00*score_all[-poscaso]

           devres[-poscaso]=sign(mres) *sqrt(-2* (mres+
			     ifelse(delta[-poscaso]==0, 0, delta[-poscaso]*log(delta[-poscaso]-mres))))

           o1<-index1[order(abs(devres[index1]),decreasing=FALSE)[1:sn1]]

           o0<-index0[order(abs(devres[index0]),decreasing=FALSE)[1:sn0]]

            pos2<-sort(c(o1,o0))


      center.trim=apply(x[pos2,],2,mean)
      sd.trim=apply(x[pos2,],2,sd)

      sd.tc<-ifelse(sd.trim==0,sd.c,sd.trim)
      center.tc<-ifelse(sd.trim==0,center,center.trim)

       x.c = sweep(x, 2, center.tc)
       xsd=sweep(x.c,2,sd.tc,"/")

        xsd.iter=xsd[pos2,]
        y.iter=Surv(t[pos2],delta[pos2])

  if(iter == "EN")
    {
     my.alpha <- seq(0.1,1,0.1)
     cvm=rep(0,length(my.alpha))
     var.selected.EN <- matrix(0,dim(xsd)[2],length(my.alpha))
     link<- matrix(0,length(pos2),length(my.alpha))

     fit.EN.cv=cva.glmnet(x=xsd.iter,y=y.iter,alpha=my.alpha,lambda=frac,family="cox",standardize=FALSE)

      for (j in 1:length(my.alpha)){

           cvm[j]=fit.EN.cv$modlist[[j]]$cvm[which(fit.EN.cv$modlist[[j]]$cvm==min(fit.EN.cv$modlist[[j]]$cvm))]

          var.selected.EN[,j]=as.matrix(coef(fit.EN.cv$modlist[[j]],s="lambda.min"))
          link[,j]=xsd.iter%*%var.selected.EN[,j]
       }


      betaff=var.selected.EN[,which(cvm==min(cvm))]
         link_EN_iter=link[,which(cvm==min(cvm))]

      } else if (iter == "LASSO") {

      stime=cv.glmnet(x=xsd.iter,y=y.iter,family='cox',alpha =0.3,standardize=FALSE)

       betaff=coef(stime,s="lambda.1se")

       link_EN_iter=as.matrix(predict(stime,newx=xsd.iter,type="link"))

        }



#logplik function
lik.cand=logplik(xsd.iter,t[pos2],delta[pos2],betaff)


pii=min(exp(Tk*(lik.cand-lik)),1)

if(rbinom(1,1,pii)==1)  {poscaso=pos2
                         likold=lik
                         lik=lik.cand
                         beta=betaff
                         link_EN=link_EN_iter


if(lik>bestlik) {
bestlik=lik
b2=beta
p2=pos2
               }       }



      }


lik=bestlik
beta=b2
pos2=p2

if(lik>max.lik)
{max.lik=lik
max.pos2=pos2
max.beta=beta
        }

result<-list(max.pos2=max.pos2, max.lik=max.lik,max.beta=max.beta)
      return(result)
       }


   statn<-as.numeric(stat)

   statmat<-matrix(0,hn,latersubset)

     for (j in 1:latersubset)
  {
   statmat[,j]=statn[betterindices[,j]]
   }


    statsum<- apply(statmat, 2, sum)
    betterindices<- betterindices[,which(statsum>2)]



   ###Cstep 2 function for 10 opitmal subsets#############

    output3=rep(0,dim(betterindices)[2])
      indices3=matrix(0,hn,dim(betterindices)[2])
      mybetas3=matrix(0,p,dim(betterindices)[2])


   for (k in 1:dim(betterindices)[2])
   {
   jieguo2<-cstep2(betterindices[,k])
   output3[k]<-jieguo2$max.lik     #likelihood function
   indices3[,k]<-jieguo2$max.pos2   #subscript of subsamples
   mybetas3[,k]<-jieguo2$max.beta # coefficient
   }

   statn<-as.numeric(stat)

   statmat<-matrix(0,hn,latersubset)

     for (j in 1:dim(indices3)[2])
  {
   statmat[,j]=statn[indices3[,j]]
   }

    statsum<- apply(statmat, 2, sum)
    indices3<- indices3[,which(statsum>2)]
    output3<- output3[which(statsum>2)]


     #keep the optimal subset with largest likelihood
    lastindices_indi=which.max(output3)
    lastindices=indices3[,lastindices_indi]
    lastcrit=output3[lastindices_indi]


    max.pos2=lastindices
    devres=rep(1,n)


      center.trim=apply(x[max.pos2,],2,mean)
      sd.trim=apply(x[max.pos2,],2,sd)
      sd.tc<-ifelse(sd.trim==0,sd.c,sd.trim)
      center.tc<-ifelse(sd.trim==0,center,center.trim)
      x.c = sweep(x, 2, center.tc)
       xsd=sweep(x.c,2,sd.tc,"/")

      xsd.raw=xsd[max.pos2,]
      y.raw=Surv(t[max.pos2],delta[max.pos2])


   if(method == "EN")
    {
       my.alpha <- seq(0.1,1,0.1)
       cvm=rep(0,length(my.alpha))

       var.selected.EN <- matrix(0,p,length(my.alpha))

      link<- matrix(0,length(max.pos2),length(my.alpha))

      fit.EN.cv=cva.glmnet(x=xsd.raw,y=y.raw,alpha=my.alpha,family="cox",standardize=FALSE)


    for (j in 1:length(my.alpha)){

           cvm[j]=fit.EN.cv$modlist[[j]]$cvm[which(fit.EN.cv$modlist[[j]]$cvm==min(fit.EN.cv$modlist[[j]]$cvm))]

          var.selected.EN[,j]=as.matrix(coef(fit.EN.cv$modlist[[j]],s="lambda.min"))
          link[,j]=xsd.raw%*%var.selected.EN[,j]

       }

       max.beta=var.selected.EN[,which(cvm==min(cvm))]
         link_EN=link[,which(cvm==min(cvm))]


      }else if (method == "LASSO") {

    stime=cv.glmnet(x=xsd.raw,y=y.raw,family='cox',standardize=FALSE)

     max.beta=coef(stime,s="lambda.1se")
     link_EN=as.matrix(predict(stime,newx=xsd.raw,type="link"))

      }

    max.lik=logplik(xsd.raw,t[max.pos2],delta[max.pos2],max.beta)
    link_raw=xsd%*%max.beta

   ### reweighted step

       strata=rep(1,length(max.pos2))

       time <-t[max.pos2]
       status <- delta[max.pos2]

	sorted <- order(strata, time)
	strata <- strata[sorted]
	newstrat <- as.integer(c(1*(diff(as.numeric(strata))!=0), 1))
      stime <- as.double(time[sorted])
      sstat <- as.integer(status[sorted])

      score=exp(link_EN)[sorted]
      weights=rep(1,length(max.pos2))
      coxres <- .C("coxmart", as.integer(length(max.pos2)),
				1,
				stime,
				sstat,
				newstrat,
				as.double(score),
				as.double(weights),
				resid=double(length(max.pos2)))


            resid <- double(length(max.pos2))
            resid[sorted] <- coxres$resid     ##Matingaile residual

          	devres[max.pos2]=sign(resid) *sqrt(-2* (resid+
			     ifelse(status==0, 0, status*log(status-resid))))

           score_all=exp(xsd%*%max.beta)

           score[score<=0.0000001]=0.0000001


           h0=(status-resid)/exp(link_EN)


           h0[h0<=0.0000001]=0.0000001

           hh<-as.data.frame(cbind(time,h0))

           myfit<-lqs(log(h0)~log(time),data=hh,method="lms")

           h00=exp(as.vector(coefficients(myfit))[1]+as.vector(coefficients(myfit))[2]*log(t[-max.pos2]))

           res=delta[-max.pos2]-h00*score_all[-max.pos2]

         devres[-max.pos2]=sign(res) *sqrt(-2* (res+
			     ifelse(delta[-max.pos2]==0, 0, delta[-max.pos2]*log(delta[-max.pos2]-res))))

        q1 <- qnorm(0.995)
        q2 <- qnorm(0.995)

       ok1<- which(devres<=q1&devres>=0)
       ok2<- which(devres>=-q2&devres<0)

       ok=c(ok1,ok2)
       rwt.indices<-sort(ok)
       outliers=(1:n)[-ok]

       rwt.indices<-sort(ok)

   statn<-as.numeric(stat)


    if (sum(statn[rwt.indices])>2)  {
     rwt.indices<-rwt.indices
      } else  {

         ind1<-which(stat==1)
         ind2<-intersect(ind1,(1:n)[-rwt.indices])

         k3<-order(abs(devres[ind2]),decreasing=FALSE)[1:3]


       rwt.indices<-c(rwt.indices,k3)
      }


       #normalize by trimmed mean and standard

       center.trim=apply(x[rwt.indices,],2,mean)
      sd.trim=apply(x[rwt.indices,],2,sd)

     sd.tc<-ifelse(sd.trim==0,sd.c,sd.trim)
     center.tc<-ifelse(sd.trim==0,center,center.trim)

     x.c = sweep(x, 2, center.tc)
     xsd=sweep(x.c,2,sd.tc,"/")
     xsd.rwt=xsd[rwt.indices,]
      y.rwt=Surv(t[rwt.indices],delta[rwt.indices])

    if(method == "EN")
    {
        my.alpha <- seq(0.1,1,0.1)
      cvm=rep(0,length(my.alpha))
      var.selected.EN <- matrix(0,p,length(my.alpha))
       link<- matrix(0,length(rwt.indices),length(my.alpha))

     fit.EN.cv=cva.glmnet(x=xsd.rwt,y=y.rwt,alpha=my.alpha,family="cox",standardize=FALSE)

      for (j in 1:length(my.alpha)){

           cvm[j]=fit.EN.cv$modlist[[j]]$cvm[which(fit.EN.cv$modlist[[j]]$cvm==min(fit.EN.cv$modlist[[j]]$cvm))]

          var.selected.EN[,j]=as.matrix(coef(fit.EN.cv$modlist[[j]],s="lambda.min"))
          link[,j]=xsd.rwt%*%var.selected.EN[,j]

       }

    rwt.betas=var.selected.EN[,which(cvm==min(cvm))]

    link_rwtb=link[,which(cvm==min(cvm))]


    }else if (method == "LASSO") {
       stime=cv.glmnet(x=xsd.rwt,y=y.rwt,family='cox',standardize=FALSE)

     rwt.betas=coef(stime,s="lambda.1se")
      link_EN=as.matrix(predict(stime,newx=xsd,type="link"))

      }

       strata=rep(1,length(rwt.indices))

       time <-t[rwt.indices]
       status <- delta[rwt.indices]

	sorted <- order(strata, time)
	strata <- strata[sorted]
	newstrat <- as.integer(c(1*(diff(as.numeric(strata))!=0), 1))

      stime <- as.double(time[sorted])
      sstat <- as.integer(status[sorted])

      score=exp(link_rwtb)[sorted]
      weights=rep(1,length(rwt.indices))
      coxres <- .C("coxmart", as.integer(length(rwt.indices)),
				1,
				stime,
				sstat,
				newstrat,
				as.double(score),
				as.double(weights),
				resid=double(length(rwt.indices)))


            resid <- double(length(rwt.indices))
            resid[sorted] <- coxres$resid     ##Matingaile residual

          	devres[rwt.indices]=sign(resid) *sqrt(-2* (resid+
			     ifelse(status==0, 0, status*log(status-resid))))

           score_all=exp(xsd%*%rwt.betas)

           score[score<=0.0000001]=0.0000001

           h0=(status-resid)/exp(link_rwtb)

           h0[h0<=0.0000001]=0.0000001

           hh<-as.data.frame(cbind(time,h0))

           myfit<-lqs(log(h0)~log(time),data=hh,method="lms")

           h00=exp(as.vector(coefficients(myfit))[1]+as.vector(coefficients(myfit))[2]*log(t[-rwt.indices]))

           res=delta[-rwt.indices]-h00*score_all[-rwt.indices]

         devres[-rwt.indices]=sign(res) *sqrt(-2* (res+
			     ifelse(delta[-rwt.indices]==0, 0, delta[-rwt.indices]*log(delta[-rwt.indices]-res))))

       q1 <- qnorm(0.995)
        q2 <- qnorm(0.995)

       ok1<- which(devres<=q1&devres>=0)
      ok2<- which(devres>=-q2&devres<0)

       ok=c(ok1,ok2)
       rwt.indices2<-sort(ok)
        outliers2=(1:n)[-ok]

    rwt.max.lik=logplik(xsd.rwt,t[rwt.indices],delta[rwt.indices],rwt.betas)


    link_rwt=xsd%*%rwt.betas



   return(list(lik=max.lik,units=max.pos2,max.beta=max.beta,link_raw=link_raw,trim=gamma,rwt.betas=rwt.betas,devres=devres,rwt.indices=rwt.indices,rwt.indices2=rwt.indices2,rwt.max.lik=rwt.max.lik,rwt.outliers=outliers,rwt.outliers2=outliers2,link_rwtb=link_rwtb,link_rwt=link_rwt))


    }
