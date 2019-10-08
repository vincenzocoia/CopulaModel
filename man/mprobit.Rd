\name{mprobit}
\Rdversion{1.1}
\alias{mprobit}
\alias{pmfmordprobit}
\alias{dvineKLfn}
\alias{dvineKLss}
\alias{ARprobitvsDvine}
\alias{rvineKLfn}
\alias{f90rvineKL}
\alias{rvineKLss}
\alias{mprobitvsRvine}
\alias{f90mprobitvsRvine}
\title{
multivariate ordinal probit and approximation by vine distributions
}
\description{
multivariate ordinal probit and approximation by vine distributions
}
\usage{
pmfmordprobit(zcuts,rmat,iprint=F,ifixseed=F) # multivariate ordinal probit pmf
dvineKLfn(parvec,ucuts,pr)  # KL divergence of multivariate probit and D-vine
rvineKLfn(parvec,ucuts,A,pr) # KL divergence of multivariate probit and R-vine
dvineKLss(parvec,ucuts,pr,iprint=F) # KL sample size of multiv probit and D-vine
rvineKLss(parvec,ucuts,A,pr,iprint=F) # KL sample size of multiv probit and R-vine
f90rvineKL(parvec,ucuts,A,M,pr) # KL divergence of multiv probit and R-vine
ARprobitvsDvine(ucuts,rmat,iprint=F,prlevel=1,mxiter=50,ifixseed=F) # best D-vine
mprobitvsRvine(ucuts,rmat,A,iprint=F,prlevel=1,mxiter=50,ifixseed=F) # best R-vine
f90mprobitvsRvine(ucuts,rmat,A,iprint=F,prlevel=1,mxiter=50,ifixseed=F)
}
\arguments{
\item{parvec}{parameter vector of partial correlations with length d*(d-1)/2}
\item{rmat}{dxd correlation matrix}
\item{ucuts}{(ncateg+1)xd matrix of cut points for ordinal, e.g.,
computed from unifcuts via ucuts=unifcuts(y), 
ucuts=rbind(rep(0,d),ucuts,rep(1,d))}
\item{zcuts}{cutpoints on N(0,1) for ordinal responses: (ncat+1)xd,
where first row is a substitute for -Inf and last row is a substitute
for +Inf}
\item{A}{dxd vine array with 1:d on diagonal}
\item{M}{dxd maximal array for vine array A}
\item{pr}{vector outputted from pmfmordprobit()}
\item{iprint}{print flag for intermediate results}
\item{mxiter}{max iterations for nlm()}
\item{ifixseed}{F by default, if T, then seed is fixed before each call to
pmvnorm within the function}
\item{prlevel}{print.level for nlm() minimization}
}
\value{vector of pmf in ordinal categories in lexicographic order, for pmfmordprobit()

Kullback-Leibler divergence for dvineKLfn and rvineKLfn

KL divergence and sample size and probability vector from vine approximation,
for dvineKLss and rvineKLsss

for ARprobitvsDvine,  
there are several components:

$mordprobitpr = probability vector from multivariate ordinal probit
(length= ncateg^d)

$dvineparam = vector of parameters for D-vine approximation with
bivariate Gaussian pair-copulas

$KLdiv = KL divergence of mult probit and D-vine approximation

$KLss = KL sample size

$vinepr = probability vector from D-vine approximation

for mprobitvsRvine and f90mprobitvsRvine, 
there are several components:

$mordprobitprmat = prob vectors from multivariate ordinal probit 
(ncateg^d x d!/2) for different permutations

$parmat = parameter of R-vine approximation for each permutation (C(d,2) x d!/2)

$vKLdiv = vector of KL divergences of multivariate probit and R-vine

$vKLss = vector of KL sample sizes

$rvineprmat= columns of probability vectors from R-vine approximation 
(ncateg^d x d!/2)
}
\details{pmfmordprobit() requires library mvtnorm for the function
pmvnorm for multivariate normal rectangle probabilities for dimensions
3 and higher. f90mprobitvsRvine is the faster version of mprobitvsRvine,
with some computations in fortran90.
f90mprobitvsRvine() and mprobitvsRvine() require library combinat.
}
\examples{
# D-vine
ncateg=3
ucuts3=matrix(c(.4,.5,.4, .7,.8,.6),ncateg-1,3,byrow=TRUE)
zcuts3=qnorm(ucuts3)
zcuts3=rbind(rep(-6,3),zcuts3,rep(6,3))
ucuts3=rbind(rep(0.00001,3),ucuts3,rep(0.99999,3))
rh=.6  
rmat3=toeplitz(rh^(0:2))
library(mvtnorm)
pmf3=pmfmordprobit(zcuts3,rmat3,iprint=TRUE)
dvineappr=nlm(dvineKLfn,c(rh,rh,0),hessian=TRUE,iterlim=40,print.level=1,ucuts=ucuts3,pr=pmf3)
dvineKLss(dvineappr$estimate,ucuts3,pmf3)
# multivariate AR probit versus discrete D-vine
ARprobitvsDvine(ucuts3,rmat3,iprint=FALSE,prlevel=1,mxiter=50)
# R-vine
d=4
A=Dvinearray(d)
out=varray2M(A); M=out$mxarray
ucuts4=matrix(c(.4,.5,.4,.3, .7,.8,.6,.6),2,4,byrow=TRUE)
ucuts4=rbind(rep(0.00001,d),ucuts4,rep(.99999,d))
parvec=c(.5,.5,.5,.1,.1,0)
pr=rep(1/81,81)
out=f90rvineKL(parvec,ucuts4,A,M,pr)
print(out)
rvineKLfn(parvec,ucuts4,A,pr)
# multivariate probit versus discrete R-vine
C3=Cvinearray(3)
ncateg=3
ucuts3=matrix(c(.4,.5,.4, .7,.8,.6),ncateg-1,3,byrow=TRUE)
ucuts3=rbind(rep(0.00001,3),ucuts3,rep(.99999,3))
rmat3=toeplitz(0.5^(0:2))
library(combinat)
out3=mprobitvsRvine(ucuts3,rmat3,C3,iprint=TRUE,prlevel=1)
out3b=f90mprobitvsRvine(ucuts3,rmat3,C3,iprint=TRUE,prlevel=1) # much faster
# these are not exactly the same output because 
# pmfmordprobit uses pmvnorm() in library mvtnorm and
#  pmvnorm() has a quasi-random component
out3=mprobitvsRvine(ucuts3,rmat3,C3,iprint=TRUE,prlevel=1,ifixseed=TRUE)
out3b=f90mprobitvsRvine(ucuts3,rmat3,C3,iprint=TRUE,prlevel=1,ifixseed=TRUE) 
# the above are the same if ifixseed=TRUE
}
\keyword{distribution}
\keyword{vine}
