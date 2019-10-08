\name{pdhessmin}
\Rdversion{1.1}
\alias{pdhessmin}
\alias{pdhessminb}
\title{
Minimization with modified Newton-Raphson and positive definite Hessian
}
\description{pdhessmin: Minimization with modified Newton-Raphson and
positive definite Hessian

pdhessminb: Minimization with modified Newton-Raphson and 
positive definite Hessian, with some parameters fixed (at bounds).
}
\usage{
pdhessmin(param,objfn,dstruct,LB,UB,mxiter=30,eps=1.e-6,bdd=5,iprint=F)
pdhessminb(param,objfn,ifixed,dstruct,LB,UB,mxiter=30,eps=1.e-6,bdd=5,iprint=F)
}
\arguments{
\item{param}{starting point for minimization of function objfn() }
\item{objfn}{objective function of form objfn(param,dstruct,iprint=F);
use iprint=T to print out extra information for debugging your function.
objfn returns a list with fnval=functionvalue, grad=gradient, hess=hessian; 
that is, objfn computes the first and second order derivatives of objfn(). }
\item{dstruct}{data structure with data sets and other variables/controls to be
passed and used by objfn() }
\item{ifixed}{logical vector of same length as param, ifixed[i]=TRUE iff 
param[i] is fixed at the given value }
\item{LB}{lower bound of components of param, usually of length(param),
could also be a scalar for a common lower bound }
\item{UB}{upper bound of components of param, usually of length(param),
could also be a scalar for a common upper bound }
\item{mxiter}{maximum number of Newton-Raphson iterations}
\item{eps}{tolerance for Newton-Raphson iterations, stop when two
consecutive iterations with eps in absolute value }
\item{bdd}{bound on difference of 2 consecutive iterations, default 5}
\item{iprint}{print flag for intermediate output for each iteration
of the Newton-Raphson method}
}
\value{
\item{parmin}{parameter value at point of minimum}
\item{fnval}{function value at the minimum}
\item{invh}{inverse Hessian at the minimum, estmated covariance matrix at
MLE if objfn is negative log-likelihood}
\item{iconv}{1 for convergence, 0 for not}
\item{iposdef}{1 for positive definite Hessian at last iteration, 0 for not}
}
\details{The algorithm is due to P Krupskii.}
\examples{
data(euro07gf)
udat=euro07gf$uscore
n=nrow(udat)
d=ncol(udat)
np=2*d
stfrk2=rep(3,np);
LB.frk2=rep(-60,np);  UB.frk2=rep(60,np);
gl=gausslegendre(15)
dstructfrk=list(copname="frank",data=udat,quad=gl,repar=0);
ifixed=rep(FALSE,np);
ml= pdhessminb(stfrk2,f90cop2nllk,ifixed=ifixed,dstruct=dstructfrk,
  LB=LB.frk2,UB=UB.frk2,iprint=TRUE,eps=1.e-4);
}
\keyword{optimize}
