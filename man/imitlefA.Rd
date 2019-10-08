\name{imitlefA}
\Rdversion{1.1}
\alias{integMittagLefflerA}
\alias{imitlefA}
\alias{pimitlefA}
\alias{dimitlefA}
\alias{pcondimitlefA}
\alias{qcondimitlefA}
\alias{pimitlefAr}
\alias{dimitlefAr}
\alias{pcondimitlefAr}
\alias{imitlefA.cpar2tau}
\title{
Bivariate Archimedean copula based on integrated Mittag-Leffler LT
}
\description{
Bivariate Archimedean copula based on integrated Mittag-Leffler Laplace
transform
}
\usage{
pimitlefA(u,v,cpar)
dimitlefA(u,v,cpar)
pcondimitlefA(v,u,cpar)  # C_{2|1}(v|u;cpar)
qcondimitlefA(p,u,cpar,eps=1e-08,mxiter=30,iprint=F) # C_{2|1}^{-1}(p|u;cpar) 
imitlefA.cpar2tau(cpar)
pimitlefAr(u,v,cpar)  # reflected/survival version of pimitlefA
dimitlefAr(u,v,cpar)
pcondimitlefAr(v,u,cpar)  
}
\arguments{
\item{u}{value in interval 0,1; could be a vector }
\item{v}{value in interval 0,1; could be a vector }
\item{p}{quantile in interval 0,1; could be a vector }
\item{cpar}{parameter vector: 2-dimensional vector (or 2-column matrix)
with parameters vth>0 and de>1.  }
\item{eps}{tolerance for convergence}
\item{mxiter}{maximum number of iterations}
\item{iprint}{print flag for iterations}
}
\value{
(conditional) cdf or pdf or quantile value(s),
or Kendall's tau.
}
\references{
Joe H (2014). Dependence Modeling with Copulas.
Chapman & Hall/CRC.
}
\examples{
u=seq(.1,.6,.1)
v=seq(.4,.9,.1)
vth=.5; de=1.6; cpar=c(vth,de)
pp=pcondimitlefA(v,u,cpar)
print(pp)
qcondimitlefA(pp,u,cpar)
cdf=pimitlefA(u,v,cpar)
pdf=dimitlefA(u,v,cpar)
tau=imitlefA.cpar2tau(cpar)
print(tau)
}
\keyword{copula}
