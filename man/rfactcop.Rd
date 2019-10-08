\name{rfactcop}
\Rdversion{1.1}
\alias{rfactcop}
\alias{r1fact}
\alias{r2fact}
\alias{rbifact}
\alias{rnestfact}
\title{
Simulation for copula models with factor structure
}
\description{
Simulation for copula models with factor structure,
with links to C code 
}
\usage{
r1fact(n,d,cpar,copcode)
r2fact(n,d,cpar,copcode)
rbifact(n,grsize,cpar,copcode)
rnestfact(n,grsize,cpar,copcode)
}
\arguments{
\item{n}{sample size}
\item{d}{number of variables}
\item{cpar}{copula parameter vector}
\item{grsize}{vector of group sizes for mgrp groups with sum(grsize)=d}
\item{copcode}{current options are 1 for Gaussian, 2 for t, 3 for Gumbel,
-3 for reflected Gumbel, 5 for Frank;
9 for BB1 1-factor, 9 for BB1/Frank 2-factor and bi-factor;
11 for Gumbel/BB1 nested factor}
}
\details{
These do simulations via C code; currently there are not as many possible 
options compared with sim1fact(), sim2fact(), simbifact() and simnestfact()
}
\value{
nxd matrix of d dependent U(0,1) variables
}
\seealso{
\code{\link{factorcopsim}}
\code{\link{structcop}}
}
\examples{ 
set.seed(123)
rhvec=c(.8,.7,.6,.5,.5)
udat1=r1fact(3,5,rhvec,1); print(udat1)
udat2=r1fact(3,5,c(rhvec,4),2); print(udat2)
cpar1=seq(1.1,1.5,.1)
udat3=r1fact(3,5,cpar1,3); print(udat3)
cparbb1=c(.2,1.1,.2,1.2,.2,1.3,.2,1.4,.2,1.5)
udat9=r1fact(3,5,cparbb1,9); print(udat9)
#
set.seed(123)
cpar2=c(seq(1.1,1.5,.1),rep(1.1,5))
udat3=r2fact(3,5,cpar2,3); print(udat3)
cpar2bb1=c(cparbb1,seq(1.1,1.5,.1))
udat9=r2fact(3,5,cpar2bb1,9); print(udat9)
#
grsize=c(2,2,2)
set.seed(123)
cparbi=c(seq(1.1,1.6,.1),rep(1.1,6))
udat3=rbifact(3,grsize,cparbi,3); print(udat3)
cparbibb1=c(cparbb1,.2,1.6,seq(1.1,1.6,.1))
udat9=rbifact(3,grsize,cparbibb1,9); print(udat9)
#
grsize=c(2,2,2)
set.seed(123)
cparne=c(rep(1.1,3),seq(1.1,1.6,.1))
udat3=rnestfact(3,grsize,cparne,3); print(udat3)
cparnebb1=c(seq(0.6,1.1,.1),cparbb1)
udat9=rnestfact(3,grsize,cparnebb1,9); print(udat9)
}
\keyword{copula}
\keyword{factor model}
\keyword{simulation}
