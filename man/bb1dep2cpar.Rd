\name{bb1dep2cpar}
\Rdversion{1.1}
\alias{bb1.dep2cpar}
\alias{bb1.tau2eqlm}
\title{
BB1 Bivariate copula: mapping of  measure of association and upper tail
dependence lambda to 2-dimensional parameter vector
}
\description{
BB1 Bivariate copula: mapping of one of (beta,tau,rhoS,rhoN) and upper tail
dependence lambda to 2-dimensional parameter vector (theta,delta);
also mapping of Kendall's tau to copula parameter vector such that lower and
upper tail dependence parameters are the same.
}
\usage{
bb1.dep2cpar(value,lmU,type="tau",iprint=F)
bb1.tau2eqlm(tau,destart=1.5,mxiter=30,eps=1.e-6,iprint=F)
}
\arguments{
\item{value}{value of dependence measure}
\item{lmU}{value of upper tail dependence lambda}
\item{type}{one of "tau", "beta", "rhoS" or "rhoN"}
\item{tau}{Kendall's tau}
\item{iprint}{print flag for intermediate results}
\item{destart}{starting point of delta for iterations, default 1.5}
\item{mxiter}{maximum number of iterations}
\item{eps}{tolerance for stopping}
}
\value{
2-dimensional copula parameter (theta, delta) if a solution is found for 
bb1.dep2cpar

vector with th, de, lm1, lmu for bb1.tau2eqlm
}
\seealso{
\code{\link{deppar2taurhobetalambda}}
}
\examples{
bb1.dep2cpar(.4,.1,"beta")
bb1.dep2cpar(.4,.1,"tau")
bb1.dep2cpar(.4,.1,"rhoS")
bb1.dep2cpar(.4,.1,"rhoN")
bb1.tau2eqlm(.4)
}

\keyword{copula}
\keyword{dependence}
