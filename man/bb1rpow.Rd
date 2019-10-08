\name{bb1rpow}
\Rdversion{1.1}
\alias{bb7pow}
\alias{pbb1rpow}
\alias{pcondbb1rpow}
\alias{qcondbb1rpow}
\alias{dbb1rpow}
\alias{pbb7pow}
\alias{pcondbb7pow}
\alias{dbb7pow}
\title{
Three-parameter bivariate copula families
}
\description{
Three-parameter bivariate copula families with a power parameter
}
\usage{
pbb1rpow(u,v,cpar)
dbb1rpow(u,v,cpar)
pcondbb1rpow(v,u,cpar)
qcondbb1rpow(p,u,cpar,pvec=c(0.01,seq(.02,.98,.02),.99),icheck=F)
pbb7pow(u,v,cpar)
dbb7pow(u,v,cpar)
pcondbb7pow(v,u,cpar)
}
\arguments{
\item{u}{value in interval 0,1; could be a vector}
\item{v}{value in interval 0,1; could be a vector}
\item{p}{quantile in interval 0,1; could be a vector}
\item{cpar}{3-parameter vector: (theta,delta) or reflected BB1 or
BB7 copula with a third parameter eta>0 for the power}
\item{pvec}{for qcondbb1rpow(), a vector of quantiles to use with monotone 
interpolation}
\item{icheck}{flag to output a check on inverse}
}
\value{
(conditional) cdf or pdf or quantile value(s)
}
\examples{
cpar=c(1.1,1.5,0.5)
cpar1=c(1.1,1.5,1)
u=seq(.1,.9,.2)
v=u
pbb1rpow(u,v,cpar)
pbb7pow(u,v,cpar)
pbb1rpow(u,v,cpar1) 
pbb1r(u,v,cpar1[1:2]) # same as above
pbb7pow(u,v,cpar1) 
pbb7(u,v,cpar1[1:2]) # same as above
}
\keyword{copula}
