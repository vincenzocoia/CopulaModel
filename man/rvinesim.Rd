\name{rvinesim}
\Rdversion{1.1}
\alias{rvinesim}
\alias{cvinesim}
\alias{dvinesim}
\alias{rvinesim0}
\alias{rvinesim1}
\alias{rvinesim2}
\alias{rvinesimvec}
\alias{rvinesimvec2}
\title{
Simulation for R-vine copulas
}
\description{Simulation for R-vine copulas,
bivariate copulas in the same parametric family for each tree;
current version assume each bivariate copula family is permutation
symmetric
}
\usage{
cvinesim(p,parmat,qcond,pcond)  # same pcond for all trees
dvinesim(p,parmat,qcond,pcond,iprint=F)  # same pcond for all trees
rvinesim0(p,A,parmat,qcond,pcond,iprint=F)  # same pcond for all trees
rvinesim1(p,A,parmat,qcondnames,pcondnames,iprint=F)
rvinesim2(p,A,parvec,np,qcondnames,pcondnames,iprint=F)
rvinesimvec(nsim,A,parvec,np,qcondnames,pcondnames,iprint=F)
  # common pair-copula family for each tree
rvinesimvec2(nsim,A,ntrunc,parvec,np,qcondmat,pcondmat,iprint=F)
  # can be different pair-copula family for each edge of vine
}
\arguments{
\item{p}{vector of length d, e.g. runif(d)}
\item{nsim}{sample size for simulation}
\item{A}{dxd vine array with 1:d on diagonal,
or ntrunc x d vine array as only ntrunc rows are used}
\item{ntrunc}{truncation level between 1 and d-1}
\item{qcond}{function for inverse conditional cdf
\eqn{C_{U|V}^{-1}(u|v)}}
\item{pcond}{function for conditional cdf \eqn{C_{U|V}(u|v)}}
\item{qcondnames}{names of inverse conditional cdf functions 
\eqn{C_{U|V}^{-1}(u|v)}, for trees 1,...,ntrunc (so ntrunc=length(qcondnames) 
for these implementations)}
\item{pcondnames}{names of conditional cdf functions \eqn{C_{U|V}(u|v)},
for trees 1,...,ntrunc}
\item{qcondmat}{matrix of names of conditional quantile functions for trees 1,...,ntrunc}
\item{pcondmat}{matrix of names of conditional cdfs for trees 1,...,ntrunc}
\item{parmat}{dxd matrix:
for rvinesim1, where all bivariate copula families have 1 parameter, 
parameter in parmat[ell,j] for ell<j is the parameter of the copula
associated with A[ell,j]}
\item{parvec}{vector: for rvinesim2 and rvinesimvec,
with the union of the parameters associated with
the copulas in A[ell,j], j=ell+1,...,d. ell=1,...,ntrunc}
\item{np}{dxd matrix: for rvinesim2 and rvinesimvec,
the dimension of the vector for 
the copulas in A[ell,j], j=ell+1,...,d. ell=1,...,ntrunc;
the function will determine parvec[ip1:ip2] for the copula
associated with A[ell,j]}
\item{iprint}{print flag for intermediate results}
}
\value{
vector of length d of values on (0,1) for cvinesim, dvinesim, rvinesim0,
rvinesim1 and rvinesim2.

data matrix (each element between 0 and 1) of dimension nsimxd for 
rvinesimvec
}
\seealso{
\code{\link{rvinenllk}}
}
\examples{ 
d=5
C=Cvinearray(d); D=Dvinearray(d)
p=c(.4,.5,.6,.4,.7)
parmat=matrix(c(0,2,3,4,5,0,0,6,6,3,0,0,0,4,3,0,0,0,3,3,0,0,0,0,2),5,5,byrow=TRUE)
cvinesim(p,parmat,qcondfrk,pcondfrk)
rvinesim0(p,C,parmat,qcondfrk,pcondfrk) # same as above line
dvinesim(p,parmat,qcondfrk,pcondfrk)
rvinesim0(p,D,parmat,qcondfrk,pcondfrk) # same as above line
#
d=7
qcondnames=rep("qcondfrk",3)
pcondnames=rep("pcondfrk",3)
parvec=c(2,3,4,5,6,6, 3,4,3,3,3, 2,2,2,2)
np=matrix(0,d,d)
np[1,2:d]=1; np[2,3:d]=1; np[3,4:d]=1
A3=vnum2array(d,300)
nsim=10
uu1=matrix(0,nsim,d)
for(i in 1:nsim)
{ pp=runif(d)
  uu1[i,]=rvinesim2(pp,A3,parvec,np,qcondnames,pcondnames,iprint=FALSE)
}
print(summary(uu1))
print(cor(uu1))
# vectorized rvinesim
uu2=matrix(0,nsim,d)
uu2=rvinesimvec(nsim,A3,parvec,np,qcondnames,pcondnames,iprint=FALSE)
print(summary(uu2))
print(cor(uu2))
# BB1/Frank
qcondnames3=c("qcondbb1","qcondfrk","qcondfrk")
pcondnames3=c("pcondbb1","pcondfrk","pcondfrk")
parvec=c(.8,2,.7,1.6,.3,2.1,.6,1.5,.9,2.3,.9,2.4, 3,4,3,3,3, 2,2,2,2)
np2=matrix(0,d,d)
np2[1,2:d]=2; np2[2,3:d]=1; np2[3,4:d]=1
uu3=rvinesimvec(nsim,A3,parvec,np2,qcondnames3,pcondnames3,iprint=FALSE)
print(uu3)
# version with different pair-copula on different edges
d=4
D=Dvinearray(d)
np=matrix(0,d,d)
np[1,2:d]=1; np[2,3:d]=1; np[3,4:d]=1
qcondnames=c("qcondgum","qcondfrk","qcondfrk")
pcondnames=c("pcondgum","pcondfrk","pcondfrk")
parvec=c(2.0,2.2,1.6, 1.5,1.5, 0.3)
set.seed(1234)
rvinesimvec(2,D,parvec,np,qcondnames,pcondnames,iprint=FALSE)
pcondmat=matrix(c("",rep("pcondgum",3),"","",rep("pcondfrk",2),
  "","","","pcondfrk",rep("",4)),4,4,byrow=TRUE)
qcondmat=matrix(c("",rep("qcondgum",3),"","",rep("qcondfrk",2),
  "","","","qcondfrk",rep("",4)),4,4,byrow=TRUE)
set.seed(1234)
rvinesimvec2(2,D,ntrunc=3,parvec,np,qcondmat,pcondmat,iprint=FALSE) #as above
}
\keyword{copula}
\keyword{simulation}
\keyword{vine}
