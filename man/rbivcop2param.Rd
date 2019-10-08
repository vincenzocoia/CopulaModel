\name{rbivcop2param}
\Rdversion{1.1}
\alias{rbivcop2param}
\alias{rbb1}
\alias{rbb2}
\alias{rbb3}
\alias{rbb4}
\alias{rbb5}
\alias{rbb6}
\alias{rbb7}
\alias{rbb8}
\alias{rbb9}
\alias{rbb10}
\title{
Bivariate 2-parameter copula families: simulation
}
\description{
Bivariate 2-parameter copula families: simulation
}
\usage{
rbb1(n,cpar,type="qcond",icheck=F)
rbb2(n,cpar,type="qcond",icheck=F)
rbb3(n,cpar,type="qcond",icheck=F)
rbb6(n,cpar,type="qcond",icheck=F)
rbb7(n,cpar,type="qcond",icheck=F)
rbb4(n,cpar,icheck=F)
rbb5(n,cpar,icheck=F)
rbb8(n,cpar,icheck=F)
rbb9(n,cpar,icheck=F)
rbb10(n,cpar,icheck=F)
}
\arguments{
\item{n}{sample size}
\item{cpar}{copula parameter: vector of length 2 }
\item{type}{"qcond" for conditional approach using \eqn{C_{2|1}^{-1}};
"mix" for stochastic representation based on mixture of max-id}
\item{icheck}{flag to output means and estimated correlation of sample}
}
\value{
nx2 matrix with values in (0,1)
}
\references{
Joe H (1997). Multivariate Models and Dependence Concepts. Chapman & Hall.
}
\seealso{
\code{\link{rLTstochrep}}
}
\examples{
n=500
set.seed(12345)
udat1=rbb1(n,c(.5,1.2),icheck=TRUE)
set.seed(12345)
udat2=rbb1(n,c(.5,1.2),type="mix",icheck=TRUE)
# another method based on Mittag-Leffler LT of the BB1 Archimedean copula
set.seed(12345)
udat3=rmbb1(n,2,c(.5,1.2))
print(summary(udat3))
print(cor(udat3)[1,2])
}

\keyword{copula}
\keyword{simulation}
