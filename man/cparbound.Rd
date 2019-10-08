\name{cparbound}
\Rdversion{1.1}
\alias{cparbound}
\alias{copparbounds}
\title{
Copula parameter bounds
}
\description{
Copula parameter bounds
}
\usage{
cparbound(copname)
data(copparbounds)
}
\arguments{
\item{copname}{string for abbreviation of copula name}
}
\value{
line number in the copparbounds object with bounds of the copula parameter(s);
the lower and upper bound(s) are printed out.
}
\details{
Functions for copula cdf, pdf, conditionals, simulation do not check on
the domain of the parameter space (or other function arguments), because
it is difficult to handle all of the possibilities on the input arguments
being a scalar, vector or matrix. These
functions can be used as arguments into functions for dependence measures
and maximum likelihood.  Use simple inputs to check how functions work.
Numerical evaluation is different from mathematics; numerical instability
(e.g., roundoffs, NaN) might occur when copula parameters are such
that dependence is very strong (nearly comonotonic).  For use of these
functions for numerical maximum likelihood, if a copula parameter bound
is Inf, use a parameter value where Kendall's tau is in the range 0.9
to 0.99.
}
\seealso{
\code{\link{pcop}}
}
\examples{
cparbound("bvncop")
cparbound("pla")
cparbound("gum")
cparbound("bb1")
cparbound("bb1rpow")
cparbound("none")
data(copparbounds)
copparbounds # to see entire list of bounds
}
