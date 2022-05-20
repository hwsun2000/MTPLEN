\name{MTPLEN}
\alias{fitter}
\title{Robust penalized Cox regression based on truncation}
\usage{
fitter(t,delta,x,gamma,poscaso=NULL,D=NULL,kmax=NULL,iter,method,rseqmax=NULL,n.multistart=NULL)
}
\arguments{
\item{t}{the time variables.}
\item{delta}{event indicators.}
\item{x}{matrix of covariates.}
\item{gamma}{trimming level.}
\item{poscaso}{optional initial solution.}
\item{D}{tuning parameters for the fitting algorithm.}
\item{kmax}{tuning parameters for the fitting algorithm.}
\item{iter}{the penalized method used in iteration.There are two selection:"EN" or "LASSO".}
\item{method}{the penalized method used for last optimal subset.There are two selection:"EN" or "LASSO".}
\item{rseqmax}{tuning parameters for the fitting algorithm.}
\item{n.multistart}{tuning parameters for the fitting algorithm.}
}
\description{
The survival data of outliers were analyzed by robust Cox regression model based on truncation, and the model was fitted by truncation of individuals with small contributions to partial likelihood, so as to obtain robust estimates unaffected by outliers.
}
