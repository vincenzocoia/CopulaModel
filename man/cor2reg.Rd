\name{cor2reg}
\Rdversion{1.1}
\alias{cor2reg}
\alias{reg2cor}
\title{
correlations to regression coefficients and vice versa for vine array 
}
\description{
correlations to regression coefficients and vice versa for vine array
(general R-vine); the matrix of regression coefficients provides a simple way
for simulation for a multivariate normal distribution with correlation
matrix with the structure of a truncated vine based on the order of
variables in the vine array.
}
\usage{
cor2reg(rr,A,iprint=F)
reg2cor(phm,A)
}
\arguments{
\item{rr}{correlation matrix, dimension d}
\item{A}{dxd vine array with 1:d on diagonal; only upper triangle is used.}
\item{phm}{matrix of regression coefficients phi[j,k], dimension d}
\item{iprint}{print flag for intermediate results}
}
\value{
\item{phm}{matrix of regression coefficients from cor2reg} 
\item{rr}{correlation matrix from reg2cor}
}
\seealso{
\code{\link{cor2pcor}}
}
\examples{
d=5
C=Cvinearray(d)
pp=matrix(0,d,d)
pp[1,2:d]=0.3
set.seed(123)
pp[1,2:d]=runif(d-1,0,1)
pp[2,3:d]=0.4
rr2=pcor2cor.cvine(pp)
print(rr2)
print(cor2reg(rr2,C,iprint=FALSE)) 
pp[3,4:d]=0.6
rr3=pcor2cor.cvine(pp)
print(rr3)
phm=cor2reg(rr3,C) 
print(phm)
print(reg2cor(phm,C))
}
\keyword{vine}
\keyword{partial correlation}
