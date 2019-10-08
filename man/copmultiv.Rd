\name{copmultiv}
\Rdversion{1.1}
\alias{copmultiv}
\alias{pmfrk}
\alias{pmgal}
\alias{pmgum}
\alias{pmhr}
\alias{mhrA}
\alias{hrcondcor}
\title{
multivariate copula cdfs 
}
\description{
multivariate copula cdfs 
}
\usage{
pmfrk(uu,cpar) # multivariate exchangeable Frank
pmgum(uu,cpar) # multivariate exchangeable Gumbel
pmgal(uu,cpar) # multivariate exchangeable Galambos
pmhr(uu,cpar) # multivariate Huesler-Reiss dimension d>=3
      # cpar is d*(d-1)/2 dimensional vector with order 12, 13, 23, 14, ...
mhrA(ww,param) # exponent function of multivariate Huesler-Reiss for d>=3
hrcondcor(thmat,j) # conditional correlation matrix in jth term of A function
}
\arguments{
\item{uu}{vector of dimension d, each element in (0,1)}
\item{cpar}{parameter of the d-variate copula}
\item{ww}{vector of dimension d, each element in (0,Inf)}
\item{param}{parameter of the multivariate distribution}
\item{thmat}{d*(d-1)/2 dimensional Huesler-Reiss parameter vector as a 
symmetric matrix}
\item{j}{integer index in 1:d, when d is nrow(thmat)}
}
\value{cdf for the pm functions; 

a non-negative value for mhrA;

a (d-1)x(d-1) correlation matrix for hrcondcor.}
\details{pmhr() requires library mvtnorm}
\seealso{
\code{\link{rectmult}}
}
\examples{
pmgum(c(.8,.9,.7),2)
pmfrk(c(.8,.9,.7),2)
pmgal(c(.8,.9,.7),3)
library(mvtnorm)
pmhr(c(.8,.9,.7),c(.5,.8,1.2))  
ww=-log(c(.8,.9,.7))
exp(-mhrA(ww,c(.5,.8,1.2))) # as above with pmhr()
thmat=corvec2mat(c(.5,.8,1.2))
hrcondcor(thmat,1)
hrcondcor(thmat,2)
hrcondcor(thmat,3)
}
\keyword{copula}
