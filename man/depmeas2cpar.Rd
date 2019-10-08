\name{depmeas2cpar}
\Rdversion{1.1}
\alias{depmeas2cpar}
\title{
Convert from a dependence measure to copula parameter
}
\description{
Convert from a dependence measure to copula parameter
}
\usage{
depmeas2cpar(values,type="beta",copname="gumbel")
}
\arguments{
\item{values}{vector of values between 0 and 1 }
\item{type}{"beta", "tau", "rhoS" or "rhoN"; see deptab }
\item{copname}{copula family: current choices are "plackett", "frank", "mtcj",
"joe", "gumbel", "galambos", "huesler-reiss". }
}
\details{
The output values are obtained by interpolation from tables.
}
\value{
copula parameters with the given dependence measure values
}

\seealso{
\code{\link{deptab}}
\code{\link{deppar2taurhobetalambda}}
\code{\link{pcinterpolate}}
}
\examples{
values=c(0.29,0.79,0.41,0.88)
depmeas2cpar(values,"beta","gumbel")
depmeas2cpar(values,"tau","galambos")
depmeas2cpar(values,"rhoS","frank")
depmeas2cpar(values,"rhoN","plackett")
}
\keyword{dependence}
