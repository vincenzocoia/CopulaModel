\name{rvinediscbvnnllk}
\Rdversion{1.1}
\alias{rvinediscbvnnllk}
\alias{rvinediscbvnfullnllk}
\title{
Negative log-likelihood for discrete R-vine with Gaussian pair-copulas
and ordinal response
}
\description{
Negative log-likelihood for discrete R-vine with Gaussian pair-copulas
and ordinal response
}
\usage{
rvinediscbvnnllk(parvec,zzdat,A)  # univariate margins converted to zzdat
rvinediscbvnfullnllk(parvec,A,xmat,yvec,nrep,ncateg) # full max likelihood
}
\arguments{
\item{parvec}{parameter vector of partial correlations with length d*(d-1)/2 }
\item{zzdat}{dimension nclx(2d) with corners of rectangle , N(0,1) scale}
\item{A}{dxd vine array with 1:d on diagonal}
\item{xmat}{nn*npred matrix, nn=nrep*ncl=d*ncl, ncl=#clusters}
\item{yvec}{integer-valued vector of length nn, values in 0:(ncateg-1) 
or 1:ncateg}
\item{nrep}{cluster size d or #repeated ordinal measures}
\item{ncateg}{number of ordinal categories}
}
\value{
negative log-likelihood
}
\seealso{
\code{\link{rvinediscrete}}
}
\examples{
data(ordinalex)
xvec=c(t(ordinalex$xx))
yvec=c(t(ordinalex$yy))
uni=ordprobit.univar(xvec,yvec,iprint=FALSE)
latentdat=mord2uu(xvec,yvec,nrep=4,uni$cutpts,uni$beta)
uudat=latentdat$uudat
zzdat=latentdat$zzdat
D4=Dvinearray(4)
param=c(.5,.5,.5,.1,.1,.1)
temz=rvinediscbvnnllk(param,zzdat,D4)
print(temz)
mlz=nlm(rvinediscbvnnllk,p=param,zzdat=zzdat,A=D4,hessian=TRUE,print.level=1)
fparam=c(uni$cutpts,uni$beta, mlz$estimate)
fnllk=rvinediscbvnfullnllk(fparam,D4,xvec,yvec,nrep=4,ncateg=3)
print(fnllk)
fnllk2=rvinediscbvnfullnllk(fparam,D4,xvec,yvec-1,nrep=4,ncateg=3)
print(fnllk2)
mlf=nlm(rvinediscbvnfullnllk,fparam,A=D4,xmat=xvec,yvec=yvec,nrep=4,ncateg=3,
   hessian=TRUE,print.level=1)
}
\keyword{maximum likelihood}
\keyword{vine}
