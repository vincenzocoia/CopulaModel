\name{makedeptable}
\Rdversion{1.1}
\alias{makedeptable}
\title{
Make table of dependence measures for a 1-parameter bivariate copula family
}
\description{
Make table of dependence measures for a 1-parameter bivariate copula family
}
\usage{
makedeptable(bvec,bfn,dcop,pcop,pcond12,pcond21,LBcpar=0,UBcpar=Inf,
   itaildep=F,lmfn,zero=0,zbd=6,iprint=F)
}
\arguments{
\item{bvec}{ vector of Blomqvist beta values in increasing order, first is 0, 
last is 1 }
\item{bfn}{function to get copula parameter cpar given Blomqvist beta}
\item{dcop}{function for copula density c}
\item{pcop}{function for copula cdf C}
\item{pcond12}{function for copula conditional \eqn{C_{1|2}(u|v)}}
\item{pcond21}{function for copula conditional \eqn{C_{2|1}(v|u)}}
\item{LBcpar}{parameter for independence}
\item{UBcpar}{parameter for comonotonicity}
\item{itaildep}{T to compute tail dependence parameter lambda=lm}
\item{lmfn}{function to compute lambda given copula parameter}
\item{zero}{0 or something like 1.e-6, used with pcond12, pcond21}
\item{zbd}{integration bound with respect to N(0,1) margins for rhoN}
\item{iprint}{print flag for intermediate results, default=F}
}
\value{
table with column headings of "cpar","beta","tau","rhoS","rhoN","lambda"
}
\seealso{
\code{\link{deptab}}
\code{\link{depmeas2cpar}}
}
\examples{
bvec=seq(0,.9,.02)  # more grid points
bvec=c(bvec,.95)
bvec=c(bvec,1)
np=length(bvec)
bvec=bvec[-(np-1)]  # 0.95 can be a problem for some copula families
np=length(bvec)
frk.deptab=makedeptable(bvec,bfn=frk.b2cpar,pcop=pfrk,
  pcond12=pcondfrk,pcond21=pcondfrk,LBcpar=0,
  UBcpar=Inf,itaildep=FALSE,zero=0,zbd=7,iprint=TRUE)
tauinteg=frk.deptab[,3]
for(i in 2:(np-1))
{ frk.deptab[i,3]=frk.cpar2tau(frk.deptab[i,1]) }
print(frk.deptab)
cat("accuracy of 2-dimensional numerical integration\n")
print(abs(tauinteg-frk.deptab[,3]))
}
\keyword{dependence}
