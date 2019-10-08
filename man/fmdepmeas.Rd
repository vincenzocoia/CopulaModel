\name{fmdepmeas}
\Rdversion{1.1}
\alias{srho1fact}
\alias{srho2fact}
\alias{rhomvn}
\alias{rho2nestfact}
\alias{cop2srho}
\title{
Spearman rho matrix for 1-factor, 2-factor, nested factor copulas
and multivariate Gaussian factor model
}
\description{
Spearman rho matrix for 1-factor, 2-factor, nested factor copulas
and multivariate Gaussian factor model
}
\usage{
srho1fact(param1,pcondcop,nq)
srho2fact(param1,param2,pcondcop1,pcondcop2,nq)
rho2nestfact(param1,param2,dcop1,pcondcop2,nq)
rhomvn(fct,param) # matrix of Spearman rho for MVN factor model
cop2srho(param,pcop,nq) # Spearman rho via Gauss-Legendre quadrature
}
\arguments{
\item{param1}{ vector/matrix of copula parameters for linking copulas to 
factor 1 or global latent variable}
\item{param2}{vector of copula parameters for linking copulas to factor 2
or group latent variable}
\item{fct}{number of factors}
\item{param}{dxfct matrix of partial correlation parameters for multivariate 
normal/Gaussian factor model with fct factors }
\item{pcondcop}{function for conditional cdf of copula family for factor 1 }
\item{pcondcop1}{function for conditional cdf of copula family for factor 1 }
\item{pcondcop2}{function for conditional cdf of copula family for factor 2
or group latent variable }
\item{dcop1}{function for bivariate copula pdf for global latent variable }
\item{pcop}{function for bivariate copula cdf}
\item{nq}{number of quadrature points for Gauss-Legendre quadrature }
}
\details{
2-dimensional Gauss-Legendre quadrature for 1-factor model,
4-dimensional Gauss-Legendre quadrature for 2-factor model.
}
\value{
dxd matrix for srho1fact, srho2fact, rhomvn

single number for rho2nestfact, cop2srho
}
\seealso{
\code{\link{factorcopsim}}
\code{\link{structcop}}
}
\examples{
n=40000  
cpar=frk.b2cpar(seq(.7,.8,.02)); cpar2=frk.b2cpar(seq(.5,.6,.02))
set.seed(123)
frk1dat=sim1fact(n,cpar,qcondfrk,"frank")
set.seed(123)
frk2dat=sim2fact(n,cpar,cpar2,qcond1=qcondfrk,qcond2=qcondfrk,"frank","frank")
# 1-factor
r1=cor(frk1dat)
print(r1)
r1t=srho1fact(cpar,pcondfrk,nq=15)
print(r1t)
print(max(abs(r1-r1t)))
# 2-factor
r2=cor(frk2dat)
print(r2)
r2t=srho2fact(cpar,cpar2,pcondfrk,pcondfrk,nq=15)
print(r2t)
print(max(abs(r2-r2t)))
# nested-factor
grsize=c(4,4)
set.seed(123)
parne=c(c(4,3),rep(6,4),rep(6.5,4))
udatne=simnestfact(n=600, grsize, cop=5, parne) # Frank pair-copulas
parwglobal=c(4,3)
parwgrp=c(6,6.5)
srh=rho2nestfact(parwglobal,parwgrp,dfrk,pcondfrk,21)
print(srh)
print(cor(udatne[,1],udatne[,5:8]))  # should be round same as above
# Spearman's rho for bivariate copula (which could be margin of factor copula)
cop2srho(c(2.5,2.8), pfact1gum, 35)  # Gauss-Legendre
rhoS(c(2.5,2.8),pfact1gum)  # cubature
cop2srho(matrix(c(2.5,2.8,1.5,1.6),2,2), pfact2gum, 35) 
}
\keyword{dependence}
\keyword{factor model}
