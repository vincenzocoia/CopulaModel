\name{contourBivCop}
\Rdversion{1.1}
\alias{contourBivCop}
\title{
Contour plot of bivariate copula density with standard Gaussian margins
}
\description{
Contour plot of bivariate copula density with standard Gaussian margins
}
\usage{
contourBivCop(cpar,zvec,dcop,irefl=F)
}
\arguments{
\item{cpar}{copula parameter}
\item{zvec}{grid points of N(0,1) density for computing bivariate density}
\item{dcop}{function for bivariate copula density with arguments
u,v,cpar}
\item{irefl}{flag to take reflection of dcop, default is F}
}
\value{
nothing is returned, a contour plot is drawn
}
\examples{
zvec=seq(-3,3,.2)
# contour density plots for Frank, Gumbel, reflected Gumbel 
contourBivCop(4,zvec,dcop=dfrk)
contourBivCop(2.5,zvec,dcop=dgum)
contourBivCop(2.5,zvec,dcop=dgum,irefl=TRUE)
}
\keyword{copula}
