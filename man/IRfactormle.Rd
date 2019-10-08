\name{IRfactormle}
\Rdversion{1.1}
\alias{IRfactormle}
\alias{ml1irfact}
\alias{ml2irfact}
\alias{f90ml1irfact}
\alias{f90ml2irfact}
\alias{f90ml2irfactb}
\alias{ir1nllk}
\alias{ir2nllk}
\alias{f90irfisherinfo1}
\alias{f90irfisherinfo2}
\title{
MLE for discrete factor models for item response
}
\description{
MLE for discrete factor copula models for item response;
f90 in function name means the log-likelihood and derivatives are computed
in fortran90 code, ir1nllk and ir2nllk are also based on f90.
}
\usage{
ml1irfact(nq,start,ydata,pcond,LB=0,UB=50,ihess=F,prlevel=0,mxiter=50)
f90ml1irfact(nq,start,ydata,copname,LB=0,UB=40,ihess=F,prlevel=0,
   mxiter=50,nu1=3)
ml2irfact(nq,start,ydata,pcond1,pcond2,LB=0,UB=40,ihess=F,prlevel=0,mxiter=50)
f90ml2irfact(nq,start,ydata,copname,LB=0,UB=40,ihess=F,prlevel=0,
   mxiter=50,nu1=5,nu2=5)
f90ml2irfactb(nq,start,ifixed,ydata,copname,LB=0,UB=40,ihess=F,prlevel=0,
   mxiter=50,nu1=5,nu2=5)
ir1nllk(param,dstruct,iprfn=F)
ir2nllk(param,dstruct,iprfn=F)
f90irfisherinfo1(param,ucutp,nq=15,copname,nu1=3,nn=1000)
f90irfisherinfo2(param,ucutp,copname,nq=15,ifixed,nu1=3,nu2=3,nn=1000)
}
\arguments{
\item{param}{parameter for ir1nllk and ir2nllk, these functions are input to 
pdhessmin or pdhessminb}
\item{nq}{number of quadrature points}
\item{start}{starting point of param for nlm optimization}
\item{ifixed}{logical vector of same length as param, ifixed[i]=TRUE iff 
param[i] is fixed at the given value}
\item{ydata}{nxd integer-valued matrix with values in 0,1,...,ncat-1,
ncat is number of ordinal categories. }
\item{pcond}{function for bivariate copula conditional cdf, linking to 
latent variable}
\item{pcond1}{function for copula conditional cdf for factor 1}
\item{pcond2}{function for copula conditional cdf for factor 2}
\item{copname}{"gumbel" or "gaussian" or "t" or "gumbelt" (Gumbel for factor
1 and t for factor 2) or "tgumbel" (t for factor 1 and Gumbel for factor 2)}
\item{dstruct}{
structure that includes $quad for the Gauss-Legendre nodes and weights,
$copname for the model, $data for ydata, $cutp for the ordinal cutpoints}
\item{LB}{lower bound of components of param, usually of length(param),
could also be a scalar for a common lower bound}
\item{UB}{upper bound of components of param, usually of length(param),
could also be a scalar for a common upper bound}
\item{ihess}{option for hessian in nlm() }
\item{prlevel}{print.level in nlm()}
\item{mxiter}{maximum number of iterations iterlim in nlm() }
\item{ucutp}{(ncat-1)xd matrix of cutpoints in (0,1)}
\item{nu1}{df parameter for factor 1 if copname="t" or "gumbelt" or "tgumbel" }
\item{nu2}{df parameter for factor 2 if copname="t" or "gumbelt" or "tgumbel"}
\item{iprfn}{flag for printing of function value and derivatives }
\item{nn}{nominal sample size for inverse of Fisher information}
}
\value{$fnval, $grad, $hess for ir1nllk and ir2nllk;
MLE etc for ml1irfact and ml2irfact.

$finfo with Fisher information matrix, $SE with sqrt(diagonal of inverse
Fisher information matrix divided by nn) for
f90irfisherinfo1 and f90irfisherinfo2.
}
\seealso{
\code{\link{IRfactorsim}}
\code{\link{factorcopmle}}
}
\references{
Nikoloulopoulos A K and Joe H (2015).
Factor copula models for item response data, Psychometrika, 80, 126-150.
}
\examples{
data(ltmconv)
d=ncol(sci)
#1-factor (3 methods)
nq=21
\dontrun{
library(abind)  # need abind() for ir1factpmf and ir2factpmf
#                 called by ml1irfact and ml2irfact respectively
ml1a=ml1irfact(nq,start=rep(2,d),sci,pcond=pcondgum,LB=1,UB=20,ihess=TRUE,
  prlevel=1,mxiter=50)
ml1b=f90ml1irfact(nq,start=rep(2,d),sci,copname="gumbel",LB=1,UB=20,ihess=TRUE,
  prlevel=1,mxiter=50)
ucutp=unifcuts(sci)
gl=gausslegendre(nq)
dstrgum=list(copname="gumbel",dat=sci,quad=gl,cutp=ucutp)  
ml1c=pdhessmin(param=rep(2,d),ir1nllk,dstruct=dstrgum,LB=rep(1,d),UB=rep(20,d),
  iprint=TRUE,eps=1.e-5);
#2-factor (2 methods)
param=c(1.5,1.1,1.6,2.5,1.05,1.2,1.5,rep(.4,d))
dfdefault=2
ml2a=ml2irfact(nq,start=param,sci,pcond1=pcondgum,pcond2=pcondt,
  LB=c(rep(1,d),rep(-1,d)),UB=c(rep(20,d),rep(1,d)),prlevel=1,mxiter=50)
dstrgumt=list(copname="gumbelt",dat=sci,quad=gl,cutp=ucutp,nu2=2) 
ml2b = pdhessmin(param,ir2nllk,dstruct=dstrgumt,
  LB=c(rep(1,d),rep(-1,d)),UB=c(rep(20,d),rep(1,d)),iprint=TRUE,eps=1.e-5); }
# Fisher information (check for near non-identifiability)
theta=c(0.5,0.6,0.5,0.6,0.4)
delta=c(0.3,0.4,0.3,0.4,0.2)
ifixed=rep(FALSE,2*d)
nq=21
K=4
d=5
ucut=seq(1:(K-1))/K
ucuts=matrix(rep(ucut,d),ncol=d)
nn=1000
ifixed=rep(FALSE,2*d)
cat("\nt5/t5\n")
finft=f90irfisherinfo2(c(theta,delta),ucutp=ucuts,copname="t",nq=nq,ifixed,
   nu1=5,nu2=5,nn=nn)
print(finft$SE)
cat("\nt20/t20\n")
finft=f90irfisherinfo2(c(theta,delta),ucutp=ucuts,copname="t",nq=nq,ifixed,
   nu1=20,nu2=20,nn=nn)
print(finft$SE)  # larger as cliser to Gaussian non-identifiable case
}
\keyword{copula}
\keyword{factor model}
\keyword{maximum likelihood}
