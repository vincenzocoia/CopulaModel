\name{rvinediscrete}
\Rdversion{1.1}
\alias{rvinediscrete}
\alias{rvinepmf.discrete}
\alias{dvinepmf.discrete}
\alias{rvinediscretenllk}
\alias{rvinediscfullnllk}
\alias{rvinepmf.ordinal}
\alias{rvineordinalnllk}
\title{
Negative log-likelihoods for regular discrete vine models
}
\description{
Probabilities and negative log-likelihoods for regular discrete vine models
}
\usage{
rvinepmf.discrete(parmat,u1vec,u2vec,A,M,pcopnames,iprint=F)
dvinepmf.discrete(parmat,u1vec,u2vec,pcopnames,iprint=F)
rvinediscretenllk(parvec,uudat,A,pcopnames,LB=0,UB=10,iprint=F)
rvinediscfullnllk(parvec,ydat,xdat,nrep,upmfcdf,npar1,A,pcopnames,
  LB=0,UB=10,zero=0.00001,one=0.99999,iprint=F) # repeated measures
rvinepmf.ordinal(parmat,yvec,A,M,pcopnames,ucuts,iprint=F)
rvineordinalnllk(parvec,ydat,ncateg,A,pcopnames,LB=0,UB=10,iprint=F)
# all of these assume 1-dimensional parameter for each edge of vine
}
\arguments{
\item{parmat}{dxd parameter matrix for the model, 1-parameter pair-copulas,
position according to A}
\item{parvec}{parameter vector for the model}
\item{u1vec}{d-dimension lower vector of hyperrectangle}
\item{u2vec}{d-dimension upper vector of hyperrectangle}
\item{uudat}{dimension nx(2d) with corners of rectangle probabilities for each discrete vector observation for rvinediscretenllk}
\item{yvec}{integer-valued d-vector, values in 0:(ncateg-1)
for rvinepmf.ordinal}
\item{ncateg}{number of categories for rvineordinalnllk}
\item{ydat}{d-dimension ordinal or discrete response; nxd matrix,
values in 0:(ncateg-1) or 1:ncateg }
\item{xdat}{covariate matrix ((n*d)xq), q=#covariates, d=#repeated measurements,
for rvinediscfullnllk}
\item{nrep}{#repeated measurements per subject for rvinediscfullnllk}
\item{ucuts}{(ncateg+1)xd matrix of cut points for ordinal, computed from unifcuts
via ucuts=unifcuts(y), ucuts=rbind(rep(0,d),ucuts,rep(1,d)), for rvinepmf.ordinal}
\item{upmfcdf}{function that computes the pmf,cdf up to the value mx}
\item{npar1}{number of parameters in upmfcdf}
\item{A}{dxd vine array with 1:d on diagonal}
\item{M}{dxd vine maximum array, $mxarray component of varray2M(A)}
\item{pcopnames}{string vector with names of pair-copula cdfs of length ntrunc,
ntrunc=truncation level}
\item{LB}{lower bound of components of param}
\item{UB}{upper bound of components of param}
\item{zero}{either 0 or epsilon depending on the copula}
\item{one}{either 1 or 1-epsilon depending on the copula}
\item{iprint}{print flag for intermediate steps}
}
\details{
dvinepmf.discrete() was written before rvinepmf.discrete() because
the algorithm looks simpler for a D-vine versus the general R-vine.
rvinepmf.discrete() can be tested with a D-vine and matched to the
output of rvinepmf.discrete().
}
\value{probability for rvinepmf.discrete and rvinepmf.ordinal

negative log-likelihood value for the nllk functions
}
\seealso{
\code{\link{discreteresponse}}
\code{\link{rvinediscbvnnllk}}
}
\references{
Panagiotelis A, Czado C and Joe H (2012).
Pair Copula Constructions for Multivariate Discrete Data.
Journal of the American Statistical Association, 107, 1063-1072.
}
\examples{
# discrete response
D=Dvinearray(4); out=varray2M(D); M=out$mxarray
parmat=matrix(0,4,4); parmat[1,2:4]=2; parmat[2,3:4]=1.7; parmat[3,4]=1.1
u1vec=c(.1,.2,.1,.2); u2vec=c(.5,.6,.5,.4)
pcopnames=rep("pgum",3)
rvinepmf.discrete(parmat,u1vec,u2vec,D,M,pcopnames,iprint=FALSE)
d=5
dd=d*(d-1)/2; n=5
uudat=matrix(c(
0.50,0.49,0.28,0.43,0.43,0.63,0.61,0.47,0.56,0.56,
0.31,0.50,0.73,0.39,0.17,0.51,0.63,0.78,0.52,0.33,
0.06,0.66,0.48,0.49,0.00,0.15,0.72,0.58,0.59,0.07,
0.00,0.00,0.00,0.00,0.00,0.58,0.57,0.55,0.52,0.53,
0.99,0.98,0.98,0.99,0.98,0.99,0.98,0.98,0.99,0.98),n,2*d,byrow=TRUE)
D5=Dvinearray(d) 
out=varray2M(D5); M=out$mxarray
rvinediscretenllk(rep(2.2,dd),uudat,D5,pcopnames=rep("pfrk",4),LB=-10,UB=30)
mle=nlm(rvinediscretenllk,p=rep(2.2,dd),hessian=TRUE,print.level=1,
  uudat=uudat,A=D5,pcopnames=rep("pfrk",4),LB=rep(-10,dd),UB=rep(20,dd))
# ordinal response
data(ltmconv)
d=ncol(sci); nitem=ncol(sci); nc=4  
perm=c(6,4,3,1,5,7,2)
sciperm=sci[,perm]
A=vnum2array(d,320); out=varray2M(A); M=out$mxarray
ucuts=unifcuts(sciperm); ucuts=rbind(rep(0,d),ucuts,rep(1,d)) 
par3=c(2,2,2,2,2,2, 1.3,1.3,1.3,1.3,1.3, 1.1,1.1,1.1,1.1)
parmat3=matrix(0,d,d)
parmat3[1,2:d]=par3[1:6]
parmat3[2,3:d]=par3[7:11]
parmat3[3,4:d]=par3[12:15]
pcopnames=rep("pgum",3)
rvinepmf.ordinal(parmat3,sciperm[16,],A,M,pcopnames,ucuts)
nllk=rvineordinalnllk(par3,sciperm,nc,A,pcopnames,LB=rep(1,15),UB=rep(10,15),iprint=FALSE)
print(nllk)
}
\keyword{maximum likelihood}
\keyword{regression}
\keyword{vine}
