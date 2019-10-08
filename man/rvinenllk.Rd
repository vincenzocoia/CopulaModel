\name{rvinenllk}
\Rdversion{1.1}
\alias{rvinenllk}
\alias{rvinenllk.trunc}
\alias{rvinenllk1.trunc}
\alias{rvinenllk.trunc2}
\alias{rvinellkv.trunc2}
\alias{rvinenllk1.nonsimpl}
\alias{rvinelogpdf}
\alias{cvinelogpdf}
\alias{dvinelogpdf}
\title{
Negative log-likelihoods for regular vine models
}
\description{
Log probability density and negative log-likelihoods for regular vine models
}
\usage{
rvinelogpdf(uvec,A,parmat,logdcop,pcond,iprint=F) # this is a template
cvinelogpdf(uvec,parmat,logdcop,pcond,iprint=F) # this is a template
dvinelogpdf(uvec,parmat,logdcop,pcond,iprint=F) # this is a template
rvinenllk.trunc(parvec,udat,A,logdcopnames,pcondnames,np,ifixed,parfixed,LB=0,UB=10)
   # common pair-copula family for each tree
#rvinellkv.trunc(parvec,udat,A,logdcopnames,pcondnames,np)
rvinenllk1.trunc(parvec,udat,A,logdcopnames,pcondnames,LB=0,UB=10) 
   # scalar parameter for each edge of vine, can be used with t copulas with
   #  fixed shape parameters
rvinenllk.trunc2(parvec,udat,A,ntrunc,logdcopmat,pcondmat,np,
  ifixed,parfixed, LB=0,UB=10)
   # can be different pair-copula family for each edge of vine
rvinellkv.trunc2(parvec,udat,A,ntrunc,logdcopmat,pcondmat,np)
rvinenllk1.nonsimpl(parvec,udat,A,logdcopnames,pcondnames,ib0fixed=F,
  b0fixed=0,iprint=F,LB=0,UB=10) 
}
\arguments{
\item{parmat}{dxd parameter matrix for the model, 1-parameter pair copulas,
position according to A}
\item{parvec}{parameter vector for the model}
\item{uvec}{d-vector of uniform scores for rvinelogpdf}
\item{udat}{nxd matrix of uniform scores for rvinenllk functions}
\item{A}{dxd vine array with 1:d on diagonal}
\item{ntrunc}{truncation level between 1 and d-1}
\item{logdcop}{function for log copula pdf; same for all edge of a vine
in rvinelogpdf}
\item{pcond}{function for copula conditional cdf; same for all edge of a vine
in rvinelogpdf}
\item{logdcopnames}{string vector with names of log copula pdfs of length 
ntrunc, ntrunc=truncation level}
\item{pcondnames}{string vector with names of copula conditional cdfs of length
ntrunc, ntrunc=truncation level}
\item{logdcopmat}{matrix of names of log copula pdfs for trees 1,...,ntrunc} 
\item{pcondmat}{matrix of names of conditional cdfs for trees 1,...,ntrunc}
\item{np}{dxd where np[ell,j] is size for parameter[ell,j]
for bivariate copula in tree ell, linking variables j and A[ell,j]}
\item{ifixed}{length equal to length(param)+length(parfixed)}
\item{parfixed}{dimension equal to sum(ifixed) for the positions where
the ifixed vector is T, parfixed[1] goes to the first fixed position}
\item{LB}{lower bound of components of parvec}
\item{UB}{upper bound of components of parvec; scalar or same length as parvec}
\item{ib0fixed}{T or F for whether intercept b0 is fixed}
\item{b0fixed}{fixed intercept when copula parameter has form
  cpar=exp(b0fixed+b1*sum(ucond)) for trees 2,3.. for rvinenllk1.nonsimpl}
\item{iprint}{print flag for intermediate calculations}
}
\value{log probability density function for rvinelogpdf, cvinelogpdf, dvinelogpdf; 

vector of log pdf or log-likelihood (one element for each row of udat, for llkv functions, to be used in Vuong's procedure 

negative log-likelihood value for the nllk functions;
}
\seealso{
\code{\link{rvinenllkderiv}}
\code{\link{rvinenllkpseud}}
\code{\link{rvinesim}}
}
\examples{
uvec=c(.1,.3,.4,.5,.7)
parmat= matrix(c(0,0,0,0,0, 1.5,0,0,0,0, 1.5,1.2,0,0,0, 1.5,1.2,1.3,0,0,
   1.5,1.2,1.3,1.4,0), 5,5)
A=vbin2array(5,6) # not C or D-vine
C=Cvinearray(5); D=Dvinearray(5)

cat("\nR-vine\n")
lpdf=rvinelogpdf(uvec,A,parmat,logdcop=logdgum,pcond=pcondgum,iprint=FALSE)
print(lpdf)

cat("\nC-vine\n")
lpdf=rvinelogpdf(uvec,C,parmat,logdcop=logdgum,pcond=pcondgum,iprint=FALSE)
print(lpdf)
lpdfc=cvinelogpdf(uvec,parmat,logdcop=logdgum,pcond=pcondgum,iprint=TRUE)

cat("\nD-vine\n")
lpdf=rvinelogpdf(uvec,D,parmat,logdcop=logdgum,pcond=pcondgum,iprint=FALSE)
print(lpdf)
lpdfd=dvinelogpdf(uvec,parmat,logdcop=logdgum,pcond=pcondgum,iprint=TRUE)

d=5
qcondnames=rep("qcondfrk",4)
pcondnames=rep("pcondfrk",4)
logdcopnames=rep("logdfrk",4)
pcondnames=rep("pcondfrk",4)
parvec=c(3.6,3.6,3.6,3.6, 1.5,1.5,1.5, 1.4,1.4, 0.3)
A=vnum2array(d,3)
set.seed(123)
nsim=20
np=matrix(0,d,d)
np[1,2:d]=1; np[2,3:d]=1; np[3,4:d]=1; np[4,5]=1
udat=rvinesimvec(nsim,A,parvec,np,qcondnames,pcondnames,iprint=FALSE)
pcondmat=matrix(c("",rep("pcondfrk",4),"","",rep("pcondfrk",3),
  "","","",rep("pcondfrk",2),"","","","","pcondfrk",rep("",5)),5,5,byrow=TRUE)
logdcopmat=matrix(c("",rep("logdfrk",4),"","",rep("logdfrk",3),
  "","","",rep("logdfrk",2),"","","","","logdfrk",rep("",5)),5,5,byrow=TRUE)
rvinenllk1.trunc(parvec,udat,A,logdcopnames,pcondnames,LB=-10,UB=30) 
rvinenllk.trunc(parvec,udat,A,logdcopnames,pcondnames,np,ifixed=rep(FALSE,10),
   parfixed=NA,LB=-10,UB=30) 
rvinenllk.trunc2(parvec,udat,A,ntrunc=4,logdcopmat,pcondmat,np,
   ifixed=rep(FALSE,10),parfixed=NA,LB=-10,UB=30) # same as above
mle=nlm(rvinenllk1.trunc,p=parvec,
  udat=udat,A=A,logdcopnames=logdcopnames,pcondnames=pcondnames,
  hessian=TRUE,iterlim=30,print.level=1,LB=-10,UB=30)
rvinellkv.trunc(mle$estimate,udat,A,logdcopnames,pcondnames,np) 
rvinellkv.trunc2(mle$estimate,udat,A,ntrunc=4,logdcopmat,pcondmat,np) 
}
\keyword{maximum likelihood}
\keyword{vine}
