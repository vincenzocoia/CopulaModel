\name{rvinenllkderiv}
\Rdversion{1.1}
\alias{rvinenllkderiv}
\alias{rvinenllkder1.trunc}
\alias{dvinenllkder1.trunc}
\alias{rvinenllkder2.trunc}
\title{
Negative log-likelihood and gradient for regular vine models
}
\description{
Negative log-likelihood and gradient for regular vine models
}
\usage{
rvinenllkder1.trunc(parvec,udat,A,logdcopdernames,pconddernames,LB=0,UB=10) 
dvinenllkder1.trunc(parvec,udat,logdcopdernames,pconddernames,LB=0,UB=10)
rvinenllkder2.trunc(parvec,udat,A,logdcopdernames,pconddernames,LB=0,UB=10) 
}
\arguments{
\item{parvec}{parameter vector for the model}
\item{udat}{nxd matrix of uniform scores for rvinenllk.trunc}
\item{A}{dxd vine array with 1:d on diagonal}
\item{logdcopdernames}{string vector with names of log copula pdfs and 
derivatives, length ntrunc, ntrunc=truncation level}
\item{pconddernames}{string vector with names of copula conditional cdfs and
derivatives, length ntrunc, ntrunc=truncation level}
\item{LB}{lower bound of components of parvec}
\item{UB}{upper bound of components of parvec; scalar or same length as parvec}
}
\details{
dvinenllkder1.trunc() was written before rvinenllkder1.trunc() because
the algorithm looks simpler for a D-vine versus the general R-vine.
rvinenllkder1.trunc() can be tested with a D-vine and matched to the
output of rvinenllkder1.trunc().
}
\value{
negative log-likelihood value with gradient as an atribute;
suitable for use with nlm.
}
\seealso{
\code{\link{rvinenllk}}
}
\examples{
d=5
A=vnum2array(d,3)
nsim=20
np=matrix(0,d,d)
# example 1
qcondnames=rep("qcondfrk",4)
pcondnames=rep("pcondfrk",4)
logdcopdernames=rep("logdfrk.deriv",4)
pconddernames=rep("pcondfrk.deriv",4)
parvec=c(3.6,3.6,3.6,3.6, 1.5,1.5,1.5, 1.4,1.4, 0.3)
set.seed(123)
np[1,2:d]=1; np[2,3:d]=1; np[3,4:d]=1; np[4,5]=1
udat=rvinesimvec(nsim,A,parvec,np,qcondnames,pcondnames,iprint=FALSE)
mle=nlm(rvinenllkder1.trunc,p=parvec,udat=udat,A=A,
  logdcopdernames=logdcopdernames,pconddernames=pconddernames,
  hessian=TRUE,iterlim=30,print.level=1,LB=-10,UB=30,check.analyticals=FALSE)
mle2=nlm(rvinenllkder1.trunc,p=parvec[1:7],udat=udat,A=A,
  logdcopdernames=logdcopdernames[1:2],pconddernames=pconddernames[1:2],
  hessian=TRUE,iterlim=30,print.level=1,LB=-10,UB=30,check.analyticals=FALSE)

# example 2
qcondnames=c("qcondbb1",rep("qcondfrk",3))
pcondnames=c("pcondbb1",rep("pcondfrk",3))
logdcopdernames=c("logdbb1.deriv",rep("logdfrk.deriv",3))
pconddernames=c("pcondbb1.deriv",rep("pcondfrk.deriv",3))
parvec=c(0.5,1.6,0.5,1.6,0.5,1.6,0.5,1.6, 1.5,1.5,1.5, 1.4,1.4, 0.3)
np[1,2:d]=2; np[2,3:d]=1; np[3,4:d]=1; np[4,5]=1
set.seed(123)
udat=rvinesimvec(nsim,A,parvec,np,qcondnames,pcondnames,iprint=FALSE)
lb=c(rep(c(0,1),4),rep(-10,6)) 
ub=c(rep(c(6,6),4),rep(30,6)) 
mle=nlm(rvinenllkder2.trunc,p=parvec,udat=udat,A=A,
 logdcopdernames=logdcopdernames,pconddernames=pconddernames,
 hessian=TRUE,iterlim=30,print.level=1,LB=lb,UB=ub,check.analyticals=FALSE)
mle2=nlm(rvinenllkder2.trunc,p=parvec[1:11],udat=udat,A=A,
 logdcopdernames=logdcopdernames[1:2],pconddernames=pconddernames[1:2],
 hessian=TRUE,iterlim=30,print.level=1,LB=lb[1:11],UB=ub[1:11],check.analyticals=FALSE)
}
\keyword{maximum likelihood}
\keyword{vine}
