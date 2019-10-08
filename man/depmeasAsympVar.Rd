\name{depmeasAsympVar}
\Rdversion{1.1}
\alias{depmeasAsympVar}
\alias{blomqvist.avar}
\alias{rhoS.avar}
\alias{ktau.avar}
\title{Asymptotic variance of
Kendall's tau, Spearman's rho and Blomqvist's beta for a bivariate copula
}
\description{Asymptotic variance of
Kendall's tau, Spearman's rho and Blomqvist's beta for a bivariate copula
}
\usage{
ktau.avar(cpar,pcop,pcond12,pcond21,zero=0,tol=1.e-5)
rhoS.avar(cpar,pcop,pcond12,pcond21,nq=25,zero=0,tol=1.e-5)
blomqvist.avar(cpar,pcop)
}
\arguments{
\item{cpar}{copula parameter: scalar or vector}
\item{pcop}{function for bivariate copula cdf}
\item{pcond12}{function for \eqn{C_{1|2}(u|v)}}
\item{pcond21}{function for \eqn{C_{2|1}(v|u)}}
\item{zero}{boundary is [zero,1-zero]^2 for integration,
default is zero, but choose something like zero=1.e-6 if there is a boundary 
problem}
\item{tol}{desired accuracy for numerical integration, default 0.00001}
\item{nq}{number of quadrature points per dimension for an inner integral.
default is 25}
}
\value{
asymptotic variance = avar, so that approximate variance for sample size
n is avar/n.
}
\seealso{
\code{\link{bivdepmeas}}
}
\examples{
cpar=2
avarkt=ktau.avar(cpar,pfrk,pcondfrk,pcondfrk)
avarsp=rhoS.avar(cpar,pfrk,pcondfrk,pcondfrk,nq=25)
avarbl=blomqvist.avar(cpar,pfrk)
print(c(cpar,avarkt,avarsp,avarbl))
}
\keyword{dependence}
