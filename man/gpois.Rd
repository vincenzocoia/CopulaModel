\name{gpois}
\Rdversion{1.1}
\alias{gpois}
\alias{gpoispmfcdf}
\alias{dgpois}
\alias{pgpois}
\title{
Generalized Poisson cdf and pmf
}
\description{
Generalized Poisson cdf and pmf
}
\usage{
gpoispmfcdf(ub,theta,vrh)
dgpois(y,param)
pgpois(y,param)
}
\arguments{
\item{param}{2-vector with theta=convolution parameter and
vrh=second parameter in (0,1)}
\item{theta}{convolution parameter >0}
\item{vrh}{second parameter between 0 and 1}
\item{y}{non-negative integer (or vector for dgpois)}
\item{ub}{non-negative integer}
}
\value{
table of size (ub+1)x3 with integer, pmf and cdf for gpoispmfcdf;

pmf for dgpois (can be vectorized);

cdf for pgpois;
}
\seealso{
\code{\link{gpoisson}}
}
\examples{
theta=2.5
vrh=0.5
gpoispmfcdf(6,theta,vrh)
dgpois(0:6,c(theta,vrh))
pgpois(6,c(theta,vrh))
}
\keyword{distribution}
