\name{gausslegendre}
\Rdversion{1.1}
\alias{gausslegendre}
\title{
Gauss-Legendre quadrature
}
\description{Gauss-Legendre quadrature nodes and weights
}
\usage{
gausslegendre(nq)
}
\arguments{
\item{nq}{number of quadrature points}
}
\value{
\item{nodes}{vector of length nq of nodes in the interval (0,1)}
\item{weights}{vector of length nq, sum of weights is 1}
}
\references{
Original source for the code is
"Stroud A and Secrest D (1966), Gaussian Quadrature Formulas, 
Prentice-Hall."
}

\examples{
out=gausslegendre(15)
# same as gauss.quad.prob(15,dist="uniform") in library(statmod)
print(sum(out$weights))  # should be 1
print(sum(out$weights*out$nodes))  # should be 0.5
}
\keyword{integration}
