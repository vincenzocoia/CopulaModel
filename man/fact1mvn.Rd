\name{fact1mvn}
\alias{fact1mvn}
\title{
Rectangle probability for multivariate normal with 1-factor structure
}
\description{
Rectangle probability for multivariate normal with 1-factor structure
}
\usage{
fact1mvn(lb,ub,load1,eps=1.e-6,B=6)
}
\arguments{
\item{lb}{vector of lower limits of integral/probability  }
\item{ub}{vector of upper limits of integral/probability }
\item{load1}{loading vector for factor 1}
\item{eps}{tolerance for numerical integration }
\item{B}{upper limit of integration and negative of lower limit}
}
\value{
rectangle probability 
}
\details{For the 1-factor dxd correlation structure, the rectangle probability
is a 1-dimensional integral. 
The positive exchangeable correlation matrix is a special case.}
\seealso{
\code{\link{exchmvn}}
}
\examples{
d=5
lb=rep(-1,d)
ub=c(2,1.5,1,1.5,2)
rho=.6
pr=exchmvn(lb,ub,rho)
cat("pr=exchmvn(lb,ub,rho)=",pr,"\n")
load1=rep(sqrt(rho),d)
pr2=fact1mvn(lb,ub,load1)
cat("pr=fact1mvn(lb,ub,load1)=",pr2,"\n")
# higher-dimensional example
load=c(.9,.8,.2,.2,.2,.2,.2,.2)
d=length(load)
uvec=seq(.3,.9,.2)
for(i in 1:length(uvec)) 
{ y=qnorm(uvec[i])
  jcdf=fact1mvn(rep(-6,d),rep(y,d),load)
  cat(y,jcdf,"\n")
}
}
\keyword{distribution}
\keyword{factor model}
