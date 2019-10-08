\name{mdiscretenllk}
\Rdversion{1.1}
\alias{mdiscretenllk}
\alias{emvndiscretenllk}
\title{
negative log-likelihood of multivariate copula with discrete margins 
}
\description{
negative log-likelihood of multivariate copula with discrete margins 
}
\usage{
mdiscretenllk(cpar,uudat,pmcop,LB,UB)
emvndiscretenllk(cpar,zzdat)  # exchangeable multivariate normal
}
\arguments{
\item{cpar}{copula parameter for pmcop}
\item{uudat}{n x (2d) matrix with uudat[,1:d] lower corner and
uudat[(d+1):(2d)] upper corner of rectangle, after fitting univariate models
and converting margin via cdf}
\item{pmcop}{function for the cdf of the d-variate copula}
\item{zzdat}{n x (2d) matrix with zzdat[,1:d] lower corner and
zzdat[(d+1):(2d)] upper corner of rectangle, after fitting univariate models
and converting margin to standard normal (for emvndiscretenllk which is
positive exchangeable multivariate normal)}
\item{LB}{lower bound vector for cpar}
\item{UB}{upper bound vector for cpar}
}
\value{negative log-likelihood}
\seealso{
\code{\link{discreteresponse}}
}
\examples{
# d=4 dimensional data with transformed univariate 
uu=matrix(c(.1,.1,.1,.1, .5,.6,.7,.8, 0,0,0,0, .4,.5,.4,.3),2,8,byrow=TRUE)
zz=qnorm(uu)
zz[zz< -6]=-6
mdiscretenllk(2,uu,pmfrk,0,30)
emvndiscretenllk(0.3,zz)
}
\keyword{maximum likelihood}
