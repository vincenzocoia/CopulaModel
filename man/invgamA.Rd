\name{invgamA}
\Rdversion{1.1}
\alias{invgamA}
\alias{pinvgamA}
\alias{dinvgamA}
\alias{logdinvgamA}
\alias{pcondinvgamA}
\alias{rminvgamA}
\alias{invgamA.cpar2tau}
\alias{invgamA.tau2cpar}
\title{
Bivariate Archimedean copula based on inverse gamma LT
}
\description{
Bivariate Archimedean copula based on inverse gamma Laplace transform
}
\usage{
pinvgamA(u,v,cpar)
dinvgamA(u,v,cpar)
pcondinvgamA(v,u,cpar)  # C_{2|1}(v|u;cpar)
rminvgamA(n,d,cpar)
logdinvgamA(u,v,cpar,pgrid=0)
invgamA.cpar2tau(cpar)
invgamA.tau2cpar(tau)
}
\arguments{
\item{u}{value in interval 0,1; could be a vector}
\item{v}{value in interval 0,1; could be a vector}
\item{cpar}{parameter: could be scalar or vector (positive-valued)}
\item{n}{sample size for ripsA, positive integer}
\item{d}{dimension}
\item{pgrid}{grid of values in (0,1) to use for monotone interpolation;
see code for the default vector when pgrid is input as 0}
\item{tau}{Kendall value in (0,1)}
}
\value{
cdf, pdf, conditional cdf, conditional quantile value(s) for pinvgamA, dinvgamA,
pcondinvgamA, qcondinvgamA respectively;

log density for logdinvgamA (use for maximum likelihood);

random d-vectors for rminvgamA;

Kendall's tau for invgamA.cpar2tau;

copula parameter for invgamA.tau2cpar
}
\references{
Joe H and Hua L (2010). Tail order and intermediate tail dependence of 
multivariate copulas.
Journal of Multivariate Analysis, v 102, 1454--1471
}
\examples{
u=seq(.1,.6,.1)
v=seq(.4,.9,.1)
cpar=.5
pp=pcondinvgamA(v,u,cpar)
print(pp)
tau=invgamA.cpar2tau(cpar)
print(tau)
set.seed(123)
udata=rminvgamA(500,d=2,cpar=2) # tau=0.5
print(taucor(udata[,1],udata[,2]))
print(semicor(udata,inscore=FALSE))
ml=nlm(bivcopnllk,p=1.5,hessian=TRUE,print.level=1,
  udat=udata,logdcop=logdinvgamA,LB=.0001,UB=10)
\dontrun{
# contour of density with N(0,1) margins
zvec=seq(-3,3,.2)
contourBivCop(2,zvec,dinvgamA)}
}

\keyword{copula}
