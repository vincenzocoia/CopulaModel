\name{mvtfact}
\Rdversion{1.1}
\alias{mvtfact}
\alias{mvtpfactnllk}
\alias{mvtbifactnllk}
\alias{mvttrifactnllk}
\alias{mvtpfact}
\alias{mvtbifact}
\alias{mvttrifact}
\title{
multivariate t with common p-factor, bi-factor and tri-factor correlation 
structures 
}
\description{
multivariate t with common p-factor, bi-factor and tri-factor correlation 
structures: negative log-likelihoods and gradients
}
\usage{
mvtpfactnllk(rhvec,tdata,df)
mvtbifactnllk(rhvec,grsize,tdata,df)
mvttrifactnllk(rhvec,grsize,sbgrsize,tdata,df)
mvtpfact(tdata,start,pfact,df,prlevel=0,mxiter=100)
mvtbifact(tdata,start,grsize,df,prlevel=0,full=T,mxiter=100)
mvttrifact(tdata,start,grsize,sbgrsize,df,prlevel=0,full=T,mxiter=150)
}
\arguments{
\item{rhvec}{vector of correlation/partial correlation parameters with 
latent variables, length is d*p for mvtpfactnllk and d*2 for mvtbifactnllk}
\item{tdata}{nxd data set, univariate margins are Student t(df)}
\item{pfact}{number of factors for mvtpfact}
\item{grsize}{vector of group sizes for the bi-factor and tri-factor models, 
length mgrp for mgrp groups with sum(grsize)=d} 
\item{sbgrsize}{vector of subgroup sizes for the tri-factor model, 
length msbgrp for msbgrp groups with sum(sbgrsize)=d; sbgrsize must be
consistent with grsize as groups are split into subgroups}
\item{df}{degree of freedom parameter for multivariate Student t;
df>300 to get multivariate Gaussian}
\item{start}{starting point for numerical maximum likelihood,
length is d*p for mvtpfact and d*2 for mvtbifact}
\item{full}{T for bi-factor, F for nested-factor as special case}
\item{prlevel}{print.level for nlm}
\item{mxiter}{max number of iterations iterlim for nlm}
}
\value{
list with
negative log-likelihood $nllk and gradient vector $lgrad
for mvtpfactnllk and mvtbifactnllk

mle object (output of nlm) for mvtpfact and mvtbifact 
}
\seealso{
\code{\link{bifct}}
\code{\link{factanal.bi}}
\code{\link{factorcopmle}}
\code{\link{structcop}}
}
\references{
Krupskii P (2014). Structured Factor Copulas and Tail Inference.
PhD thesis, University of British Columbia.
}
\details{When parameters are  converted to a dxp matrix (p=2 for
bi-factor and p=3 for tri-factor),
the correlations or partial correlations in any column are unique
up to sign. 

There is further non-uniqueness for the p-factor
model with p>=2, as the loading matrix can be rotated.
To get a parameter matrix with one 0 in the second column, two 0s in the
third column etc., use pcor2load() to convert to a loading matrix,
then something like grotate2() and grotate3() to rotate the
loading matrix, then convert back to a partial correlation representation
with load2pcor().

For bi-factor and tri-factor, there is non-uniquess in the second (or
third) column if some group (subgroup) sizes are 1 or 2.

However, if the numerial optimization converges from different starting
points, the final negative log-likelihood should be the same even if
the point of convergence is not and the gradient is not close to a zero
vector.
}
\examples{ 
data(euro07gf)
udat=euro07gf$uscore
d=ncol(udat)
st1=rep(0.4,d)
st2=rep(0.4,2*d)
for(df in c(5,15))
{ tdata=qt(udat,df)
  cat("\ndf=", df,"\n")
  cat("1-factor MVt\n")
  out1t=mvtpfact(tdata,st1,pfact=1,df=df,prlevel=1)
  cat("\n2-factor MVt\n")
  out2t=mvtpfact(tdata,st2,pfact=2,df=df,prlevel=1)
  st1=out1t$estimate
  st2=out2t$estimate
}
# non-uniqueness for 2-factor
st2=matrix(st2,ncol=2)
load2=pcor2load(st2)
load2rot=grotate2(load2,row=1)
st2b=c(rep(.7,d),rep(.2,d))
out2b=mvtpfact(tdata,st2b,pfact=2,df=15,prlevel=1)
load2b=pcor2load(matrix(out2b$estimate,ncol=2))
load2brot=grotate2(load2b,row=1)
print(max(abs(load2b-load2)))
print(max(abs(load2brot-load2rot)))
print(out2t$min-out2b$min)
# bi-factor and nested-factor
df=10
tdata=qt(udat,df)
grsize=c(4,3)
bif=mvtbifact(tdata, c(rep(.8,d),.2,.2,.9,.2,.2,.8,.2),grsize,df=df,
  prlevel=1,full=TRUE,mxiter=100)
nestf=mvtbifact(tdata, c(.9,.2,rep(.8,d)),grsize,df=df,
  prlevel=1,full=FALSE,mxiter=100)
# tri-factor, simulated example
grsize=c(6,6)
sbgrsize=c(3,3,3,3)
d=sum(grsize)
p=3
tripar=((d*p):1)/(d*p+1)
param1=tripar[1:d]
param2=tripar[(d+1):(2*d)]
param3=tripar[(2*d+1):(3*d)]
n=100
df=10
robj=trifct(grsize,sbgrsize,param1,param2,param3)
achol=chol(robj$fctmat)
set.seed(123)
z=matrix(rnorm(n*d),n,d)
z=z%*%achol
udata=uscore(z)
tdata=qt(udata,df)
out=mvttrifactnllk(tripar,grsize,sbgrsize,tdata,df)
print(out)
ml=mvttrifact(tdata,start=tripar,grsize,sbgrsize,df,prlevel=1,mxiter=150)
st2=tripar; tripar[1:d]=.8
ml2=mvttrifact(tdata,start=st2,grsize,sbgrsize,df,prlevel=1,mxiter=150)
}
\keyword{factor model}
\keyword{maximum likelihood}
