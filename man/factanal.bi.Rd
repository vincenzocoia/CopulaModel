\name{factanal.bi}
\Rdversion{1.1}
\alias{pfactnllk}
\alias{bifactnllk}
\alias{trifactnllk}
\alias{factanal.co}
\alias{factanal.bi}
\alias{factanal.tri}
\title{
Gaussian factor analysis: common and bi-factor
}
\description{
Gaussian factor analysis: common and bi-factor
}
\usage{
pfactnllk(rhvec,Robs,nsize)
bifactnllk(rhvec,grsize,Robs,nsize)
trifactnllk(rhvec,grsize,sbgrsize,Robs,nsize)
factanal.co(factors,start,data=1,cormat=NULL,n=100,prlevel=0,mxiter=100)
factanal.bi(grsize,start,data=1,cormat=NULL,n=100,prlevel=0,mxiter=100)
factanal.tri(grsize,sbgrsize,start,data=1,cormat=NULL,n=100,prlevel=0,mxiter=150)
}
\arguments{
\item{factors}{the number of factors p for factanal.co}
\item{rhvec}{For pfactnllk: vector of length d*p of partial correlations, 
first d are correlations with factor 1, then
partial correlations with factor k given previous factors for k=2,...,p.
For bifactnllk: vector of length d*2; first d correlations with common factor, 
then partial correlations with group factor given common factor.
For trifactnllk: vector of length d*3; first d correlations with common factor, 
then partial correlations with group factor given common factor, and
finally partial correlations with subfactor given common and groups factors.}
\item{Robs}{dxd correlation matrix}
\item{nsize}{sample size}
\item{grsize}{vector of group sizes for the bi-factor and tri-factor models}
\item{sbgrsize}{vector of subgroup sizes for the tri-factor model}
\item{start}{starting point for nlm() to find p-factor or bi-factor
or tri-factor parameters to minimize the Gaussian negative log-likelihood}
\item{data}{nsize x d data set to compute the correlation matrix}
\item{cormat}{dxd correlation matrix}
\item{n}{sample size}
\item{prlevel}{print.level for nlm}
\item{mxiter}{maximum number of iterations for nlm}
}
\value{
For pfactnllk, bifactnllk and trifactnllk, a value for negative Gaussian
log-likelihood and the gradient vector with respect to rhvec.

For factanal.co, a list with
$nllk, $rhmat= dxp matrix of partial correlations, 
$loading=dxp loading matrix after varimax,
$rotmat=pxp rotation matrix used by varimax.
Note that the loading matrix and matrix of partial correlations
are not unique for p>=2.

For factanal.bi, a list with
$nllk, $rhmat= dx2 matrix of correlations and partial correlations.

For factanal.tri, a list with
$nllk, $rhmat= dx3 matrix of correlations and partial correlations.
}
\details{The output of factanal.co should be similar to that of
factanal(). The algorithm is different so that it may fail to converge
depending on the starting point (when p is larger).

pfactnllk, bifactnllk and trifactnllk return the gradient vector in order that
gradient calculations are faster and maybe
the number of iterations in nlm() can be reduced .
}
\seealso{
\code{\link{bifct}}
\code{\link{mvtfact}}
}
\examples{
rhpar=c( 0.81,  0.84,0.84, 0.54,0.57,0.49, 0.51,0.54,0.55,0.70,
  0.53,0.56,0.53,0.67,0.70)
Robs=corvec2mat(rhpar)
n=100
d=6
rhvec=c(rep(.5,d),rep(.4,d))
pfactnllk(rhvec,Robs,n)
factanal.co(2,start=rhvec,cormat=Robs,n=100,prlevel=1)
# bi-factor
grsize=c(3,3)
bifactnllk(rhvec,grsize,Robs,n)
factanal.bi(grsize,start=rhvec,cormat=Robs,n=100,prlevel=1)
# tri-factor
grsize=c(4,4)
sbgrsize=c(2,2,2,2)
d=sum(grsize)
p=3
tripar=((d*p):1)/(d*p+1)
param1=tripar[1:d]
param2=tripar[(d+1):(2*d)]
param3=tripar[(2*d+1):(3*d)]
n=50
robj=trifct(grsize,sbgrsize,param1,param2,param3)
achol=chol(robj$fctmat)
set.seed(123)
z=matrix(rnorm(n*d),n,d)
z=z%*%achol
zdata=nscore(z,iopt=TRUE)
Robs=cor(zdata)
trifactnllk(tripar,grsize,sbgrsize,Robs,n)
factanal.tri(grsize,sbgrsize,start=tripar,cormat=Robs,n=50,prlevel=1,mxiter=150)
}
\keyword{factor model}
