\name{bivcopnllk}
\Rdversion{1.1}
\alias{bivcopnllk}
\alias{bivcopnllk.ipol}
\alias{bivmodnllk}
\title{
Negative log-likelihood for bivariate copula model
}
\description{
Negative log-likelihood for bivariate copula model to input to 
numerical minimizer
}
\usage{
bivcopnllk(param,udat,logdcop,ivect=T,LB=0,UB=1000)
bivmodnllk(param,xdat,logdcop,logpdf1,cdf1,np1,logpdf2,cdf2,np2,ivect=T,LB,UB)
bivcopnllk.ipol(param,udat,logbpdf,logupdf,uquant,iunivar,ppvec=NULL,LB,UB)
  # latter function used interpolation for a copula based on
  #  F_{12}(F^{-1}(u1),F^{-1}(u2)) where F^{-1} is computationally difficult
}
\arguments{
\item{param}{parameter of model: copula parameter only for bivcopnllk;
(par1,par2,par3) for bivmodnllk where par1 is the vector of parameters for
univariate margin 1, par2 is for univariate margin 2, and
par3 is for the bivariate copula. }
\item{udat}{nx2 matrix, assume each column is U(0,1) distributed}
\item{xdat}{nx2 matrix, original untransformed data}
\item{logdcop}{function with log of copula density}
\item{ivect}{flag that is T if logdcop can take vectorized inputs}
\item{logpdf1}{function with log of univariate density of first variable}
\item{cdf1}{function with univariate cdf of first variable}
\item{np1}{dimension of par1 (parameters for first variable)}
\item{logpdf2}{function with log of univariate density of second variable}
\item{cdf2}{function with univariate cdf of second variable}
\item{np2}{dimension of par2 (parameters for second variable)}
\item{logbpdf}{function with log of density of \eqn{F_{12}}}
\item{logupdf}{function with log of univariate marginal density}
\item{uquant}{function for inverse of univariate marginal cdf}
\item{iunivar}{vector with indices such that par1=param[iunivar] is 
the univariate parameter}
\item{ppvec}{vector of quantiles to use for interpolation, a default is used
in this is input as NULL; a possibility is something like
ppvec=seq(min(udat),max(udat),length=100)}
\item{LB}{vector of lower bounds on parameter values}
\item{UB}{vector of upper bounds on parameter values}
}
\value{
negative log-likelihood 
}
\seealso{
\code{\link{ordinal.bivcop}}
}
\examples{
data(alae600)
alae=alae600$alae 
loss=alae600$loss
ppareto=function(x,param)
{ alp=param[1]; s=param[2]
  u=1-(1+(x/s))^{-alp}
  u
}
logdpareto=function(x,param)
{ alp=param[1]; s=param[2]
  lpdf=log(alp/s)-(alp+1)*log(1+(x/s))
  lpdf
}
paretonllk=function(param,xdat)
{ alp=param[1]; s=param[2]
  if(alp<=0 | s<=0) return(1.e10)
  xs=xdat/s
  nllk=(alp+1)*log(1+xs) - log(alp/s)
  sum(nllk)
}
alae.pareto=nlm(paretonllk,p=c(2.39,1.6),hessian=TRUE,print.level=1,xdat=alae/10000)
loss.pareto=nlm(paretonllk,p=c(1.5,1.5),hessian=TRUE,print.level=1,xdat=loss/10000)
mle1p=alae.pareto$estimate
mle2p=loss.pareto$estimate
ualae.pareto=ppareto(alae/10000,mle1p)
uloss.pareto=ppareto(loss/10000,mle2p)
udat=cbind(ualae.pareto,uloss.pareto)
# IFM or inference functions for margins
ifm.pp=nlm(bivcopnllk,p=1.5,hessian=TRUE,print.level=1,
  udat=udat,logdcop=logdgum,LB=1,UB=20)
# mle is 1.434, nllk is -80.41
# full likelihood
parampp=c(mle1p,mle2p,ifm.pp$estimate)
full.pp=nlm(bivmodnllk,p=parampp,hessian=TRUE,print.level=1,
    xdat=cbind(alae/10000,loss/10000),logdcop=logdgum,
    logpdf1=logdpareto,cdf1=ppareto,np1=2,logpdf2=logdpareto,cdf2=ppareto,np2=2,
    LB=c(rep(0,4),1),UB=c(rep(100,4),20) )
# Pareto margin with bivariate Gaussian copula as a comparison 
ml=nlm(bivcopnllk,p=.5,hessian=TRUE,print.level=1,
  udat=udat,logdcop=logdbvncop,LB=-1,UB=1)
# mle is 0.471, nllk is -72.25
logdnorm=function(x,par1) { dnorm(x,log=TRUE) }
uqnorm=function(p,par1) { qnorm(p) }
ppvec=seq(min(udat),max(udat),length=100)
mla=nlm(bivcopnllk.ipol,p=.5,udat,hessian=TRUE,print.level=1,
  logbpdf=logdbvn,logupdf=logdnorm,uquant=uqnorm,iunivar=NULL,
  ppvec=ppvec,LB=-1,UB=1)
# mle is 0.469, nllk is -72.09 with the interpolation
}
\keyword{copula}
\keyword{maximum likelihood}
