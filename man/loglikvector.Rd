\name{loglikvector}
\Rdversion{1.1}
\alias{vuongllkr}
\alias{vuong2llkr}
\alias{mdiscretellkv}
\alias{emvndiscretellkv}
\alias{rvinediscretellkv}
\alias{ir1factpmf}
\alias{ir2factpmf}
\alias{strfactllkv}
\alias{mvtbifactllkv}
\alias{mvttrifactllkv}
\alias{rvinellkv.trunc}
\title{
Vector of log-likelihoods for a model 
}
\description{
Vector of log-likelihoods for a model for use with Vuong's procedure
}
\usage{
vuongllkr(llkv1,llkv2)
vuong2llkr(llkv1,llkv2,dim1,dim2)
mdiscretellkv(param,uudat,mrectpr)  # multivariate discrete
emvndiscretellkv(param,zzdat)  # exchangeable discretized multivariate normal 
rvinediscretellkv(param,uudat,A,pcopnames,iprint=FALSE) # discrete R-vine
ir1factpmf(param,dstruct,pcondcop)
ir2factpmf(param,dstruct,pcondcop1,pcondcop2)
strfactllkv(param,udat,strmodel,copname,nq,grsize=0,nu=0,ipdf=1)
mvtbifactllkv(param,tdata,grsize,df,full=TRUE)
mvttrifactllkv(param,tdata,grsize,sbgrsize,df,full=TRUE)
rvinellkv.trunc(parvec,udat,A,logdcopnames,pcondnames,np)
}
\arguments{
\item{llkv1}{vectors of log-likelihoods for model 1, length n is same
as that for the sample size}
\item{llkv2}{vectors of log-likelihoods for model 2, length n is same
as that for the sample size}
\item{dim1}{parameter vector dimension for model 1}
\item{dim2}{parameter vector dimension for model 2}
\item{param}{parameter vector for the model}
\item{parvec}{parameter vector for the model}
\item{uudat}{dimension nx(2d) with corners of rectangle probabilities for each 
discrete vector observation on U(0,1) scale; uudat[,1:d]<uudat[(d+1):(2*d)]}
\item{zzdat}{dimension nx(2d) with corners of rectangle probabilities for each 
discrete vector observation on N(0,1) scale; zzdat[,1:d]<zzdat[(d+1):(2*d)]}
\item{mrectpr}{function for multivariate rectangle probability}
\item{A}{dxd vine array with 1:d on diagonal}
\item{pcopnames}{string vector with names of copula cdfs of length ntrunc, 
ntrunc=truncation level}
\item{pcondcop}{pcond function 1-factor ordinal model}
\item{pcondcop1}{pcond function first factor of 2-factor ordinal model}
\item{pcondcop2}{pcond function second factor of 2-factor ordinal model}
\item{dstruct}{structure with $dat for dataset, $cutp for cutpounts on (0,1)
scale, and $quad for Gauss-Legendre object with quadrature points and nodes}
\item{udat}{nxd matrix of uniform scores}
\item{strmodel}{one of "1factor", "2factor", "bifactor", "nestedfactor"}
\item{copname}{something like "frank", "bb1", "bb1frank", "bb1frk", "t", 
"tapprox"}
\item{nq}{number of quadrature points }
\item{grsize}{vector of group sizes for mgrp groups with sum(grsize)=d}
\item{sbgrsize}{vector of subgroup sizes by partitioning grsize vector}
\item{nu}{Student t df parameters of copname is "t"}
\item{ipdf}{default to 1 to compute log-likelihood only with gradients}
\item{tdata}{nxd matrix of Student t scores}
\item{df}{df parameter for multivariate Student t}
\item{iprint}{print flag for intermediate calculations}
\item{full}{T for bi-factor and F for nested-factor structure for multivariate t}
\item{logdcopnames}{string vector of names of logcopula densities for trees 1,...,ntrunc}
\item{pcondnames}{string vector of names of cond cdfs for trees 1,...,ntrunc}
\item{np}{dxd where np[ell,j] is size for parameter th[ell,j]
for bivariate copula in tree ell, variables j and A[ell,j]}
}
\value{
vector of log-likelihood values at parameter estimate, one for each 
observation, for the llkv functions; 

vector of discrete probabilities or likelihoods 
for ir1factpmf and ir2factpmf;

95 percent interval for the mean of llkv2-llkv1 for the vuongllkr and
vuong2llkr function (the latter adjusts for the Schwarz/BIC correction);
negative interval means that model 1 is better, 
positive interval means that model 2 is better, and
interval that includes 0 implies models not significantly different.
}
\seealso{
\code{\link{bivcopnllk}}
\code{\link{factorcopmle}}
\code{\link{IRfactormle}}
\code{\link{mdiscretenllk}}
\code{\link{mvtfact}}
\code{\link{rvinediscrete}}
\code{\link{structcop}}
}
\examples{
# example with discrete ordinal data
data(ltmconv)
d=ncol(sci)
n=nrow(sci)
ucutp=unifcuts(sci)
nq=21
gl = gausslegendre(nq)
# factor models
library(abind)
dstrsci=list(dat=sci,quad=gl,cutp=ucutp)
par1.gum=c(1.467402,1.071322,1.589827,2.476544,1.049539,1.191183,1.506965)
llkv1fgum=ir1factpmf(par1.gum,dstrsci,pcondgum)
llkv1fgum=log(llkv1fgum)
par2.gumt=c(1.3632887,1.5646628,1.1770374,1.3801555,1.5733537,
 1.7975962,1.2583718,
-0.3123371,0.4745622,-0.5367812,-0.6948827,0.5530323,
 0.3288671,-0.4507251)
dfdefault=2
llkv2fgumt=ir2factpmf(par2.gumt,dstrsci,pcondgum,pcondt)
llkv2fgumt=log(llkv2fgumt)
# truncated vine model
perm=c(6,4,3,1,5,7,2)  
sciperm2=sci[,perm]
d=ncol(sciperm2)
n=nrow(sciperm2)
A2=vnum2array(d,320)
out=varray2M(A2) 
M2=out$mxarray
ucuts2=unifcuts(sciperm2)
ucuts2=rbind(rep(0,d),ucuts2,rep(1,d))
# gumbel/t(5)
pbvtcop1=function(u,v,rh,df=dfdefault)
{ param=c(rh,df)
  u[u>=1]=.9999999; v[v>=1]=.9999999
  u[u<=0]=.0000001; v[v<=0]=.0000001
  xt=qt(u,df); yt=qt(v,df); 
  pbvt(xt,yt,param)
}
pcopnames2=c("pgum","pbvtcop1")
dfdefault=5
par2.rvine=c(1.1317132,1.0337752,1.2786283,1.4919012,1.3498679,1.5146284,
 0.4852748,0.1611857,-0.1258033,0.3406340,0.2977515)
parmat2=matrix(0,d,d)
parmat2[1,2:d]=par2.rvine[1:6]
parmat2[2,3:d]=par2.rvine[7:11]
llkvrvine2t=rep(0,n)
for(i in 1:n) 
{ llkvrvine2t[i]=rvinepmf.ordinal(parmat2,sciperm2[i,],A2,M2,pcopnames2,ucuts2) 
}
llkvrvine2t=log(llkvrvine2t)
cmpa=vuongllkr(llkvrvine2t,llkv1fgum)
print(cmpa)
cmpb=vuongllkr(llkvrvine2t,llkv2fgumt)
print(cmpb)

# example with continuous data
bevec=c(.8,.7,.6,.5,.5)
cpar.frk=frk.b2cpar(bevec)
cpar.gum=gum.b2cpar(bevec)
set.seed(123)
udat=sim1fact(100,cpar.frk,qcondfrk,"frank")
out.frk=ml1fact(nq=21,cpar.frk,udat,dfrk,LB=-30,UB=30,prlevel=0,mxiter=50)
out.gum=ml1fact(nq=21,cpar.gum,udat,dgum,LB=1,UB=30,prlevel=0,mxiter=50)
llkvfrk=strfactllkv(out.frk$estimate,udat,"1factor","frank",nq=21,ipdf=1)
llkvgum=strfactllkv(out.gum$estimate,udat,"1factor","gumbel",nq=21,ipdf=1)
print(c(out.frk$min,-sum(llkvfrk)))
print(c(out.gum$min,-sum(llkvgum)))
cmp=vuongllkr(llkvgum,llkvfrk)
print(cmp) # second model is "true" model so interval should be positive,
}
\keyword{maximum likelihood}
