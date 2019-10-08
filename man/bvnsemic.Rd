\name{bvnsemic}
\Rdversion{1.1}
\alias{bvnsemic}
\title{
Semi-correlation for bivariate Gaussian
}
\description{
Semi-correlation for bivariate Gaussian
}
\usage{
bvnsemic(rho)
}
\arguments{
\item{rho}{(vector of) correlation parameter for bivariate Gaussian}
}
\value{
semi-correlation (same for lower and upper)
}
\seealso{
\code{\link{rhoNsemic}}
}
\examples{
rho=seq(.1,.9,.1)
bvnsemic(rho)  
# see function for the formula
}
\keyword{dependence}
