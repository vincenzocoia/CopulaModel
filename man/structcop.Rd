\name{structcop}
\Rdversion{1.1}
\alias{structcop}
\alias{simnestfact}
\alias{simbifact}
\alias{f90str1nllk}
\alias{f90str2nllk}
\title{
Simulation and maximum likelihood for structured factor copulas
}
\description{
Simulation and maximum likelihood for structured nested and bi-factor copulas
}
\usage{
simnestfact(nn,grsize,cop,param)
simbifact(nn,grsize,cop,param)
f90str1nllk(param,dstruct,iprfn=F) # nested-factor
f90str2nllk(param,dstruct,iprfn=F) # bi-factor
}
\arguments{
\item{nn}{sample size }
\item{grsize}{vector of group sizes for mgrp groups with sum(grsize)=d} 
\item{cop}{number code for a copula model:
1 for Gaussian/normal, 2 for Student t, 3 for Gumbel, 5 for Frank;
10 for Gumbel/BB1 for simnestfact;
9 for BB1/Frank for simbifact}
\item{param}{parameter vector; length is d+mgrp+1(cop==2)}
\item{dstruct}{structure that includes $quad for the gauss-legendre nodes and 
weights, $copname for the model, $data for data set of dependent U(0,1),
$grsize for grsize. 
For t-factor model, also $nu for degree of freedom parameter.
Also $repar is a code for reparametrization (check examples).}
\item{iprfn}{flag for printing of function value and derivatives }
}
\value{data matrix of dimension nn x d for simnestfact() and simbifact()
with U(0,1) or N(0,1) or t(df) margins.

$fnval, $grad, $hess for f90str1nllk and f90str2nllk.
}
\seealso{
\code{\link{factorcopmle}}
\code{\link{factorcopsim}}
\code{\link{mvtfact}}
}
\references{
Krupskii P and Joe H (2015).
Structured factor copula models: theory, inference and computation.
J Multivariate Analysis, 138, 19-33.
}
\details{
f90str1nllk: "t", "tbb1" (t for group latent to global latent,
BB1 for observed to group latent), "tgum", "frank", "gumbel", "frkgum", 
"frkbb1", "gumbb1", "tgum", 

f90str2nllk: "frank", "gumbel", "gumfrk", "bb1frk", "bb1gum", "t"
"tapprox" (latter uses monotone interpolation for the Student t cdf)

The order of BB1 parameters for all of the models with a BB1 component
is theta1,theta2,...thetad,delta1,....deltad,
with thetas>0 and deltas>1.

This is different for BB1 parameters for the 1-factor BB1 and 2-factor
BB1/Frank, where it is theta1,delta1,theta2,delta2,...,thetad.deltad.

The difference is due to how to handle the f90 code.
}
\examples{ 
gl=gausslegendre(25)
grsize=c(4,4,3)
d=sum(grsize)
n=500
# nested-factor copula
mgrp=length(grsize)
set.seed(123)
parne=c(rep(4,3),rep(6,4),rep(6.5,4),rep(7,3))
udatne=simnestfact(n,grsize,cop=5,parne)
dstrfrk=list(data=udatne,copname="frank",quad=gl,repar=0,grsize=grsize)
npar=mgrp+d
outn= pdhessminb(rep(3,npar),f90str1nllk, ifixed=rep(FALSE,npar), dstrfrk, 
  LB=rep(0,npar), UB=rep(30,npar), mxiter=30, eps=5.e-5,iprint=TRUE)
# bi-factor copula
set.seed(123)
parbi=c(rep(4,11),rep(6,4),rep(6.5,4),rep(7,3))
udatbi=simbifact(n,grsize,cop=5,parbi)
npar=2*d
dstrfrk=list(data=udatbi,copname="frank",quad=gl,repar=0,grsize=grsize,pdf=0)
nllk=f90str2nllk(parbi,dstrfrk)
outb=pdhessminb(c(rep(2,d),rep(3,d)),f90str2nllk,ifixed=rep(FALSE,npar),dstrfrk,
  LB=rep(0,npar), UB=rep(20,npar), mxiter=30, eps=5.e-5,iprint=TRUE)
}
\keyword{copula}
\keyword{factor model}
\keyword{maximum likelihood}
\keyword{simulation}
