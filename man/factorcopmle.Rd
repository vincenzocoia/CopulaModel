\name{factorcopmle}
\Rdversion{1.1}
\alias{factorcopmle}
\alias{ml1fact}
\alias{ml1factb}
\alias{ml2fact}
\alias{ml2factb}
\alias{f90ml1fact}
\alias{f90ml2fact}
\alias{f90ml2factb}
\alias{f90cop1nllk}
\alias{f90cop2nllk}
\title{
MLE and negative log-likelihood for factor copula models
}
\description{
MLE and negative log-likelihood for factor copula models (both
common and structured factors).
f90 in function name means the log-likelihood and derivatives are computed
in fortran90 code; note that ml2fact is same as f90ml2fact.
}
\usage{
ml1fact(nq,start,udata,dcop,LB=0,UB=1.e2,prlevel=0,mxiter=100)
ml1factb(nq,start,ifixed,udata,dcop,LB=0,UB=1.e2,prlevel=0,mxiter=100)
f90ml1fact(nq,start,udata,copname,LB=0,UB=40,ihess=F,prlevel=0,mxiter=100,nu=3)
f90ml2fact(nq,start,udata,copname,LB=0,UB=40,repar=0,ihess=F,prlevel=0,mxiter=100,nu1=3,nu2=3)
f90ml2factb(nq,start,ifixed,udata,copname,LB=0,UB=40,ihess=F,prlevel=0,mxiter=100,nu1=3,nu2=3)
f90cop1nllk(param,dstruct, iprfn=F) # 1-factor
f90cop2nllk(param,dstruct, iprfn=F) # 2-factor
}
\arguments{
\item{param}{
parameter for f90cop1nllk and f90cop2nllk, these functions are input to 
pdhessmin or pdhessminb}
\item{nq}{number of quadrature points }
\item{start}{starting point of param for nlm optimization }
\item{ifixed}{logical vector of same length as param, ifixed[i]=TRUE iff 
param[i] is fixed at the given value }
\item{udata}{nxd matrix of values in (0,1) }
\item{dcop}{function for bivariate copula density }
\item{copname}{"frank", "gumbel", "gumfrank" (see below) }
\item{dstruct}{structure that includes $quad for the gauss-legendre nodes and 
weights, $copname for the model, and $data of form udata 
For t-factor model, also $nu for degree of freedom parameter.
Also $repar as a code for reparametrization (check examples).
}
\item{LB}{lower bound of components of param, usually of length(param),
could also be a scalar for a common lower bound }
\item{UB}{upper bound of components of param, usually of length(param),
could also be a scalar for a common upper bound }
\item{nu}{degree of freedom parameter for 1-factor model with t copulas}
\item{nu1}{degree of freedom parameter for first factor of 2-factor model with t copulas}
\item{nu2}{degree of freedom parameter for second factor 2-factor model with t copulas}
\item{repar}{repar=1 for reparametrization for Gumbel copula (cpar->1+param^2)
and repar=2 for BB1 copula with second parameter (delta->1+param^2)} 
\item{ihess}{option for hessian in nlm() }
\item{prlevel}{print.level in nlm()}
\item{mxiter}{max number of iterations or iterlim in nlm() }
\item{iprfn}{flag for printing of function value and derivatives }
}
\value{$fnval, $grad, $hess for f90cop1nllk and f90cop2nllk;

MLE etc for ml1fact and ml2fact.
}
\seealso{
\code{\link{factorcopsim}}
\code{\link{mvtfact}}
\code{\link{structcop}}
}
\references{
Krupskii P and Joe H (2013).
Factor copula models for multivariate data.
Journal of Multivariate Analysis, 120, 85-101.
}
\details{
ml1fact (ml1factb) uses only R code;
f90ml1fact, f90ml2fact, f90ml2factb link to Fortran 90 code;
ml2fact (ml2factb) is the same as f90ml2fact (f90ml2factb).

Current options are the following (but user can use the f90 code as
templates for adding other linking copulas).

(1a) f90ml1fact: "frank", "gumbel", "t", "bb1"

(1b) f90cop1nllk: "frank", "lfrank", "gumbel", "t", "bb1"

(2) f90cop2nllk and f90ml2fact (f90ml2factb): 
"frank", "lfrank", "gumbel", "gumfrank", "bb1frank", "t",
"tapprox" (latter uses monotone interpolation for the Student t cdf). 

For BB1, the order of parameters is theta[1],delta[1],...,theta[d],delta[d],
where thetas>0 and deltas>0.
}
\examples{
# 1-factor
cpar.frk=c(12.2,3.45,4.47,4.47,5.82); d=5
n=300
set.seed(123)
frkdat=sim1fact(n,cpar.frk,qcondfrk,"frk")
cat("\nFrank 1-factor MLE: standalone R and then f90\n")
out.frk=ml1fact(nq=21,cpar.frk,frkdat,dfrk,LB=-30,UB=30,prlevel=1,mxiter=100)
gl21=gausslegendre(21)
dstrfrk=list(copname="frk",data=frkdat,quad=gl21,repar=0)
out=pdhessminb(cpar.frk,f90cop1nllk,ifixed=rep(FALSE,d),dstruct=dstrfrk,
  LB=rep(-30,d),UB=rep(30,d),iprint=TRUE,eps=1.e-5)
set.seed(123)
dfdefault=5
tdat=sim1fact(n,c(.6,.7,.6,.8,.65),qcondt,"t")
dstrt=list(copname="t",data=tdat,quad=gl21,repar=0,nu=10)
out=pdhessminb(rep(.6,d),f90cop1nllk,ifixed=rep(FALSE,d),dstruct=dstrt,
  LB=rep(-1,d),UB=rep(1,d),iprint=TRUE,eps=1.e-5)
dstrt=list(copname="t",data=tdat,quad=gl21,repar=0,nu=5)
out=pdhessminb(rep(.6,d),f90cop1nllk,ifixed=rep(FALSE,d),dstruct=dstrt,
  LB=rep(-1,d),UB=rep(1,d),iprint=TRUE,eps=1.e-5)
# 2-factor
d=7; be1=c(.7,.6,.7,.6,.7,.6,.7); be2=c(.4,.4,.4,.4,.3,.3,.3)
cpar1.frk=frk.b2cpar(be1); cpar2.frk=frk.b2cpar(be2)
n=300
set.seed(123)
frkdat=sim2fact(n,cpar1.frk,cpar2.frk,qcondfrk,qcondfrk,"frk","frk")
cat("\nFrank 2-factor MLE: nlm and then pdhessmin\n")
out.frk=ml2fact(nq=21,c(cpar1.frk,cpar2.frk),frkdat,copname="frank",
  LB=-30,UB=30,prlevel=1,mxiter=100)
dstrfrk=list(copname="frank",data=frkdat,quad=gl21,repar=0)
out=pdhessminb(c(cpar1.frk,cpar2.frk),f90cop2nllk,ifixed=rep(FALSE,2*d),
  dstruct=dstrfrk, LB=rep(-20,2*d),UB=rep(20,2*d),iprint=TRUE,eps=1.e-5)
}
\keyword{copula}
\keyword{factor model}
\keyword{maximum likelihood}
