\name{bivdepmeas}
\Rdversion{1.1}
\alias{bivdepmeas}
\alias{ktau}
\alias{rhoS}
\alias{rhoN}
\title{
Kendall's tau, Spearman's rho and normal scores correlation for a bivariate
copula
}
\description{
Kendall's tau, Spearman's rho and normal scores correlation for a bivariate
copula
}
\usage{
ktau(cpar,icond=T,pcond12,pcond21,zero=0,dcop,pcop,B=6,tol=0.0001)
rhoS(cpar,cop,zero=0,icond=F,tol=0.0001)
rhoN(cpar,icond=T,pcond,icdf=F,pcop,dcop,B=6,tol=0.0001)
}
\arguments{
\item{cpar}{copula parameter: scalar or vector}
\item{icond}{if T, numerical integration via conditional cdf of copula;
otherwise via the copula cdf for rhoS and copula cdf/pdf for tau, rhoN}
\item{icdf}{if T, numerical integration via the cdf and Hoeffding's identity,
otherwise via the pdf (for rhoN)}
\item{pcond12}{function for \eqn{C_{1|2}(u|v)}}
\item{pcond21}{function for \eqn{C_{2|1}(v|u)}}
\item{pcond}{function for \eqn{C_{2|1}(v|u)}, for rhoN}
\item{cop}{for rhoS, copula conditional cdf if icond=T, and bivariate copula 
cdf if icond=F}
\item{zero}{boundary is \eqn{[zero,1-zero]^2} for integration,
default is 0, but choose smething like zero=1.e-6 if there is a boundary 
problem}
\item{dcop}{function for bivariate copula density, needed if icond=F}
\item{pcop}{function for bivariate copula cdf, needed if icond=F}
\item{B}{upper limit for integration over \eqn{[0,B]^2} with respect to normal
density if icond=F}
\item{tol}{desired accuracy for numerical integration, default 0.0001}
}
\value{
\item{depmeas}{dependence measure with value in (-1,1)}
}
\seealso{
\code{\link{depmeasAsympVar}}
}
\examples{
cpar=2
tau=ktau(cpar,pcond12=pcondpla,pcond21=pcondpla,zero=0)
spear=rhoS(cpar,cop=pfrk,zero=0)
rhoz=rhoN(cpar,pcond=pcondfrk)
}
\keyword{dependence}
