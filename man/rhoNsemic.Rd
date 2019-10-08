\name{rhoNsemic}
\Rdversion{1.1}
\alias{rhoNsemic}
\title{
Semi-correlations for bivariate copula
}
\description{
Semi-correlations for bivariate copula 
}
\usage{
rhoNsemic(cpar,dcop,pcop,pcond,B=6,isym=F,iinfbd=T,iprint=F)
}
\arguments{
\item{cpar}{copula parameter: scalar or vector}
\item{dcop}{function for bivariate copula pdf}
\item{pcop}{function for bivariate copula cdf}
\item{pcond}{function for \eqn{C_{2|1}(v|u)}}
\item{B}{bound on numerical integration on \eqn{[0,B]^2} with respect to
standard normal densities}
\item{isym}{flag for reflection symmetric, default is F}
\item{iinfbd}{flag for infinity as limit; if T upper bound is Inf for 
1-dimensional integral, otherwise B; default is T}
\item{iprint}{print flag for intermediate results}
}
\value{vector of length 2 with
lcorr = lower normal scores semi-correlation of bivariate copula,
ucorr = upper normal scores semi-correlation of bivariate copula.
}
\seealso{
\code{\link{pcop}}
\code{\link{semicor}}
}
\examples{
rhoNsemic(2,dfrk,pfrk,pcondfrk,B=6,isym=TRUE)
rhoNsemic(2,dgum,pgum,pcondgum,B=6,isym=FALSE)
rhoNsemic(2,dgumr,pgumr,pcondgumr,B=6,isym=FALSE)
}
\keyword{dependence}
