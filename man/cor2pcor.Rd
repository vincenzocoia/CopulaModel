\name{cor2pcor}
\Rdversion{1.1}
\alias{cor2pcor}
\alias{pcor2cor}
\alias{cor2pcor.cvine}
\alias{cor2pcor.dvine}
\alias{cor2pcor.rvine}
\alias{pcor2cor.cvine}
\alias{pcor2cor.dvine}
\alias{pcor2cor.rvine}
\alias{pcor2cor.1tr}
\alias{pcor2cor.2tr}
\alias{pcor2cor.3tr}
\alias{pcor2cor.truncvine}
\alias{vineResidVar}
\title{
correlations to partial correlations and vice versa for vine array 
}
\description{
correlations to partial correlations and vice versa for vine array
(C vine, D vine or general R vine) 
}
\usage{
cor2pcor.cvine(rr)
pcor2cor.cvine(pc)
cor2pcor.dvine(rr)
pcor2cor.dvine(pc)
cor2pcor.rvine(rr,A)
pcor2cor.rvine(pc,A,byrow=T)
pcor2cor.1tr(pp,A)  # 1-truncated vine
pcor2cor.2tr(pp,A)  # 2-truncated vine
pcor2cor.3tr(pp,A,iprint=F)  # 3-truncated vine
pcor2cor.truncvine(pp,A,ntrunc,iprint=F)
vineResidVar(rmatobj,A,ntrunc)
}
\arguments{
\item{rr}{correlation matrix, dimension d}
\item{pc}{
partial correlation array, dimension d, stored in C-vine or D-vine format}
\item{pp}{
partial correlation array, dimension d, stored in format of vine array A}
\item{byrow}{TRUE if partial correlations by tree are stored by rows in pc,
FALSE if partial correlation in location (j,k) of pc is \eqn{\rho_{j,k;S}} for
a conditioning set S.}
\item{A}{dxd vine array with 1:d on diagonal; only upper triangle is used}
\item{ntrunc}{truncation level of vine; only first ntrunc rows of A are used}
\item{rmatobj}{object from  pcor2cor.truncvine() with $rmat and $phmat}
\item{iprint}{print flag for intermediate results}
}
\value{
\item{pc}{from cor2pcor.cvine and cor2pcor.dvine, partial correlation array}
\item{pcobj}{from cor2pcor.rvine, list with $pctree and $pcmat,
the former with partial correlations by tree in rows, the latter
with \eqn{\rho_{j,k;S}} in location (j,k).
}
\item{rr}{correlation matrix from pcor2cor.cvine, pcor2cor.dvine, pcor2cor.rvine}
\item{rmatobj}{component $rmat for correlation matrix and $phmat for
matrix of regression coefficients from pcor2cor.truncvine}
\item{psi2}{vector of residual variances from vineResidVar}
}
\details{
pcor2cor.1tr, pcor2cor.2tr, pcor2cor.3tr are specific versions
of pcor2cor.truncvine; pcor2cor.1tr is called by pcor2cor.truncvine
and the other two are absorbed into pcor2cor.truncvine.
}
\seealso{
\code{\link{cor2reg}}
}
\examples{
d=5
D=Dvinearray(d)
C=Cvinearray(d)
A=vnum2array(d,bnum=3)  
# bnum in 0 to 2^dpow-1 where dpow=2^dcase, dcase=(d-2)*(d-3)/2
rr=toeplitz(c(1,.5,.25,.125,.05))
cor2pcor.dvine(rr)
cor2pcor.rvine(rr,D)
cor2pcor.cvine(rr)
cor2pcor.rvine(rr,C)
pcobj=cor2pcor.rvine(rr,A)
pcor2cor.rvine(pcobj$pctree,A)
pcor2cor.rvine(pcobj$pcmat,A,byrow=FALSE)
rmatobj=pcor2cor.truncvine(pcobj$pctree,A,ntrunc=3)
print(rmatobj$phmat)
print(rmatobj$rmat)
vineResidVar(rmatobj,A,3)
}
\keyword{vine}
\keyword{partial correlation}
