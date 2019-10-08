\name{coptrivmxid}
\Rdversion{1.1}
\alias{coptrivmxid}
\alias{pmxid3ps}
\alias{pmxidr3ps}
\alias{pmxid2ps}
\alias{pmxidr2ps}
\alias{pmxid3ls}
\alias{pmxid2ls}
\alias{ptmm1}
\alias{pbmm1}
\alias{pcondbmm1}
\alias{dbmm1}
\title{
Trivariate copula cdfs from mixture of max-id
}
\description{
Trivariate copula cdfs from mixture of max-id and bivariate (1,2) margin.
These are functions for copula-based Markov order2 time series.
}
\usage{
pmxid3ps(uu,param)   # positive stable LT and pfrk for H12,H32
pmxidr3ps(uu,param)  # reflected pmxid3ps
pmxid3ls(uu,param)   # log series LT and pfrk for H12,H32
ptmm1(uu,param)      # positive stable LT and pgum for H12,H32 (example of MM1)
pmxid2ps(u,v,param)  # bivariate (1,2) margin
pmxidr2ps(u,v,param)
pmxid2ls(u,v,param)
pbmm1(u,v,param)
pcondbmm1(v,u,param) # conditional cdf of pbmm1 with respect to u
dbmm1(u,v,param)     # density of pbmm1
}
\arguments{
\item{uu}{vector of size 3 with numbers in (0,1) for trivariate copula}
\item{u}{value in interval 0,1 for bivariate copula}
\item{v}{value in interval 0,1 for bivariate copula}
\item{param}{copula parameter =(LTpar,cpar12)}
}
\value{
cdf or pdf
}
\examples{
u=seq(.1,.9,.2)
v=u
param=c(1.5,1.2)
pmxid2ps(u,v,param)  
pmxidr2ps(u,v,param)  
pmxid2ls(u,v,param)  
pbmm1(u,v,param)  
pcondbmm1(v,u,param)  
uu=c(.1,.4,.6)
pmxid3ps(uu,param)  
pmxidr3ps(uu,param)  
pmxid3ls(uu,param)  
ptmm1(uu,param)  
}
\keyword{copula}
