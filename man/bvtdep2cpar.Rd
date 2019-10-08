\name{bvtdep2cpar}
\Rdversion{1.1}
\alias{bvt.dep2cpar}
\title{
Bivariate t copula: mapping of  measure of association 
to copula parameter given shape parameter nu>0
}
\description{
Bivariate t copula: mapping of  measure of association 
to copula parameter given shape parameter nu>0
}
\usage{
bvt.dep2cpar(value,nu,type="tau",ngrid=15,iprint=F)
}
\arguments{
\item{value}{value of dependence measure}
\item{nu}{positive shape parameter for t distribution}
\item{type}{one of "tau", "beta","lambda", "rhoS" or "rhoN"}
\item{ngrid}{number of grid points to use in interpolation, default 15}
\item{iprint}{print flag for intermediate results}
}
\value{parameter in interval (-1.1)
}
\seealso{
\code{\link{deppar2taurhobetalambda}}
}
\examples{
bvt.dep2cpar(.4,5,"beta")
bvt.dep2cpar(.4,5,"tau")
bvt.dep2cpar(.4,5,"lambda")
bvt.dep2cpar(.4,5,"rhoS")
bvt.dep2cpar(.4,5,"rhoN")
}

\keyword{copula}
\keyword{dependence}
