\name{ipsA}
\Rdversion{1.1}
\alias{integpstableA}
\alias{pipsA}
\alias{dipsA}
\alias{logdipsA}
\alias{pcondipsA}
\alias{qcondipsA}
\alias{pipsAr}
\alias{dipsAr}
\alias{logdipsAr}
\alias{pcondipsAr}
\alias{qcondipsAr}
\alias{ripsA}
\alias{ipsA.cpar2tau}
\alias{ipsA.tau2cpar}
\title{
Bivariate Archimedean copula based on integrated positive stable LT
}
\description{
Bivariate Archimedean copula based on integrated positive stable Laplace
transform
}
\usage{
pipsA(u,v,cpar)
dipsA(u,v,cpar)
pcondipsA(v,u,cpar)  # C_{2|1}(v|u;cpar)
qcondipsA(p,u,cpar)  # C_{2|1}^{-1}(p|u;cpar)
ripsA(n,cpar)
logdipsA(u,v,cpar)
ipsA.cpar2tau(cpar)
ipsA.tau2cpar(tau,mxiter=20,eps=1e-06, cparstart=0, iprint=F)
}
\details{
For the reflected copula, the functions are
pipsAr, dipsAr, pcondipsAr, qcondipsAr, logdipsAr.
}
\arguments{
\item{u}{value in interval 0,1; could be a vector}
\item{v}{value in interval 0,1; could be a vector}
\item{p}{quantile in interval 0,1; could be a vector}
\item{cpar}{parameter: could be scalar or vector (positive-valued)}
\item{n}{sample size for ripsA, positive integer}
\item{tau}{tau in interval (-1,1), could be a vector}
\item{mxiter}{maximum number of Newton-Raphson iterations}
\item{eps}{tolerance for convergence of Newton-Raphson iterations}
\item{cparstart}{starting point for Newton-Raphson iterations}
\item{iprint}{print flag for Newton-Raphson iterations}
}
\value{
cdf, pdf, conditional cdf, conditional quantile value(s) for pipsA, dipsA,
pcondipsA, qcondipsA respectively;

log density for logdipsA;

random pairs for ripsA;

Kendall's tau for psA.cpar2tau;

copula parameter for ipsA.tau2cpar
or parameter value for a given Kendall's tau.
}
\references{
Joe H and Ma C (2000). Multivariate survival functions with a min-stable
property. Journal of Multivariate Analysis, 75, 13-35.
}
\examples{
u=seq(.1,.6,.1)
v=seq(.4,.9,.1)
cpar=.6
pp=pcondipsA(v,u,cpar)
vv=qcondipsA(pp,u,cpar)
print(pp)
print(vv)
tau=ipsA.cpar2tau(cpar)
print(tau)
tauv=seq(-.9,.9,.1)
cpar=ipsA.tau2cpar(tauv)
print(cpar)
set.seed(123)
udata=ripsA(500,cpar=cpar[15]) # tau=0.5
print(taucor(udata[,1],udata[,2]))
print(cor(udata,method="kendall"))
}

\keyword{copula}
