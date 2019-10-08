\name{rcopulaGARCH}
\Rdversion{1.1}
\alias{rcopulaGARCH}
\alias{rgarch1fact}
\alias{rgarch2fact}
\alias{rgarchbifact}
\alias{rgarchnestfact}
\alias{rgarchbifactmvt}
\alias{rgarchnestfactmvt}
\title{
Simulation for copula GARCH models with factor structure
}
\description{
Simulation for copula GARCH models with factor structure,
with links to C code
}
\usage{
rgarch1fact(n,garchpar,cpar,sigma0,copcode)
rgarch2fact(n,garchpar,cpar,sigma0,copcode)
rgarchbifact(n,grsize,garchpar,cpar,sigma0,copcode)
rgarchnestfact(n,grsize,garchpar,cpar,sigma0,copcode)
rgarchbifactmvt(n,grsize,garchpar,cpar,sigma0,copcode=2) # monotone interpolation
rgarchnestfactmvt(n,grsize,garchpar,cpar,sigma0,copcode=2) # monotone interpolation
}
\arguments{
\item{n}{sample size}
\item{garchpar}{6xd matrix, where d is the number of assets,
rows are mu,ar1,omega,alpha1,beta1,nu; set ar1 to be 0 vector if AR term not
used in GARCH model}
\item{cpar}{copula parameter vector}
\item{grsize}{vector of group sizes for mgrp groups with sum(grsize)=d}
\item{sigma0}{d-vector with starting values for conditional SDs (from GARCH
output)}
\item{copcode}{current options are 1 for Gaussian, 2 for t, 3 for Gumbel,
-3 for reflected Gumbel, 5 for Frank;
9 for BB1 1-factor, 9 for BB1/Frank 2-factor and b-factor;
11 for Gumbel/BB1 nested factor} 
}
\details{
rgarchbifactmvt and rgarchnestfactmvt use monotone interpolation for the
Student t univariate cdf, for faster computations
}
\value{
\item{lgret}{dxn matrix of log returns}
\item{portfret}{nx1 portfolio return vector, assuming equally weighted}
}
\seealso{
\code{\link{rfactcop}}
\code{\link{structcop}}
}
\examples{ 
garchpar=matrix(c(0.094, 0.125, 0.068, 0.212, 0.092,
 -0.062,-0.040,-0.082,-0.020,-0.051,
  0.019, 0.014, 0.016, 0.082, 0.014,
  0.071, 0.075, 0.088, 0.092, 0.080,
  0.909, 0.914, 0.882, 0.842, 0.902,
 9.4, 8.3,9.5, 7.9,9.6), 6,5,byrow=TRUE)
sigma0=c(1.113,1.286,0.857,1.176,0.968)
garchpar6=cbind(garchpar,c(0.1, -0.05,0.01,0.08,0.9,10.))
sigma6=c(sigma0,1)
grsize=c(2,2,2)
cpar1=seq(1.1,1.5,.1)
cpar2=c(seq(1.1,1.5,.1),rep(1.1,5))
cparbi=c(seq(1.1,1.6,.1),rep(1.1,6))
cparne=c(rep(1.1,3),seq(1.1,1.6,.1))
#
set.seed(123)
out=rgarch1fact(3,garchpar,cpar1,sigma0,copcode=3)
print(out)
set.seed(123)
out=rgarch2fact(3,garchpar,cpar2,sigma0,copcode=3)
print(out)
set.seed(123)
out=rgarchbifact(3,grsize,garchpar6,cparbi,sigma6,copcode=3)
print(out)
set.seed(123)
out=rgarchnestfact(3,grsize,garchpar6,cparne,sigma6,copcode=3)
print(out)
}
\keyword{copula}
\keyword{simulation}
