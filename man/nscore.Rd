\name{nscore}
\Rdversion{1.1}
\alias{nscore}
\alias{nscoreOpta}
\title{
Transform to normal scores
}
\description{
Transform each variable in a data matrix to normal scores
and finding the adjustment so that mean of squares of scores is 1.
}
\usage{
nscore(data,iopt=F)
nscoreOpta(n,mxiter=20, eps=1.e-4,iprint=F)
}
\arguments{
\item{data}{data matrix or data frame }
\item{iopt}{
iopt=T to use adjustment 'a' in nscoreOpta
that makes normal scores have mean of 0 and
variance of 1; default is iopt=F in which case a=-0.5 is used }
\item{n}{sample size}
\item{mxiter}{maximum number of iterations}
\item{eps}{tolerance for stopping}
\item{iprint}{print flag for iterations}
}
\value{
data matrix with same number of columns as input for nscore

number near -0.5 for nscoreOpta
}
\seealso{
\code{\link{uscore}}
}
\examples{
set.seed(123)
x=matrix(rnorm(40),20,2)
z=nscore(x)
nscoreOpta(100) # -0.5766102
nscoreOpta(1000) # -0.5522775
}
\keyword{nonparametric}
