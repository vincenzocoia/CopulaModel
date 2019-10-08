\name{dbvn}
\Rdversion{1.1}
\alias{dbvn}
\alias{dbvn2}
\alias{logdbvn}
\alias{dbvt}
\title{
bivariate normal (Gaussian) and t densities
}
\description{
bivariate normal (Gaussian) and t densities, and logarithm of densities
}
\usage{
dbvn(x,rho)
dbvn2(x1,x2,rho)
logdbvn(x1,x2,rho)
dbvt(x1,x2,param,log=FALSE)
}
\arguments{
\item{x}{2-dimensional vector}
\item{x1}{scalar or vector}
\item{x2}{scalar or vector}
\item{rho}{correlation in (-1,1)}
\item{param}{(rho,nu) with rho in (-1,1) and nu>0}
\item{log}{T to return log density, =F to return density}
}
\value{
bivariate normal density (could be vector for dbvn2 if one or more of
x1,x2,rho is a vector)
}
\examples{
dbvn(c(1,.5),.5)
dbvn2(seq(.2,1,.2),seq(.1,.5,.1),.5)
dbvt(1,.5,c(.5,50))
jpdf=dbvt(1,.5,c(.5,5))
dbvt(1,.5,c(.5,5),log=TRUE)
log(jpdf)
}
\keyword{distribution}
