\name{uscore}
\Rdversion{1.1}
\alias{uscore}
\title{
Transform to U(0,1) scores
}
\description{
Transform each variable in a data matrix to uniform scores
}
\usage{
uscore(data,aunif=-0.5)
}
\arguments{
\item{data}{data matrix or data frame }
\item{aunif}{rank i -> (i+aunif)/(n+1+2*aunif) for sample size n;
default is -0.5}
}
\value{
data matrix with same number of columns as input
}
\seealso{
\code{\link{nscore}}
}
\examples{
set.seed(123)
xx=matrix(rnorm(40),20,2)
uu=uscore(xx)
}
\keyword{nonparametric}
