\name{ordinal}
\Rdversion{1.1}
\alias{ordinal}
\alias{d2v}
\alias{v2d}
\alias{d2b}
\alias{unifcuts}
\alias{ordinal2fr}
\alias{bprobitnllk}
\alias{bprobitwPrednllk}
\alias{polychoric.bivtab}
\alias{polychoric0}
\alias{polychoric}
\alias{polychoric.wPred}
\title{
multivariate ordinal response
}
\description{
multivariate ordinal response with no covariates: some utility functions
}
\usage{
d2v(d,ncat,ii,izero=F)  # decimal to category vector
v2d(jj,ncat)  # vector to decimal
d2b(d,ii)  # decimal to binary vector
unifcuts(odat)
ordinal2fr(odat,ncat)
bprobitnllk(rho,zcuts,bfr,jj1,jj2)
bprobitwPrednllk(rho,zzdat,jj1,jj2)
polychoric.bivtab(bivtab,iprint=F,prlevel=0)
polychoric0(odat,iprint=F,prlevel=0)   # number of categories can vary
polychoric(odatfr,zcuts,iprint=F,prlevel=0) # same categories for all variables
polychoric.wPred(zzdat,iprint=F,prlevel=0)
}
\arguments{
\item{d}{dimension d of ordinal response vector}
\item{ncat}{#categories (assumed labeled 0,1,...,(ncat-1) or 1,...,ncat)}
\item{rho}{latent correlation parameter in (-1,1)}
\item{ii}{non-negative integer in 0 to ncat^2-1}
\item{jj}{d-vector, each element in 0,1,...,(ncat-1)}
\item{izero}{if T, categories start at 0, otherwise 1}
\item{odat}{nxd matrix of ordinal responses in 0,...,(ncat-1) or 1,...,ncat}
\item{odatfr}{nx(d+1) matrix: d columns of ordinal responses and final
column with frequency of each distinct observed d-vector}
\item{zcuts}{cutpoints on N(0,1) for ordinal responses: (ncat+1)xd,
where first row is a substitute for -Inf and last row is a substitute
for +Inf}
\item{bivtab}{bivariate table of counts for 2 ordinal variables}
\item{zzdat}{nx(2d) matrix with corners of rectangle for each vector
observation in N(0,1) scale; for ordinal data, 
ordprobit.univar and mord2uu can been used to get zzdat}
\item{bfr}{vector of bivariate frequencies}
\item{jj1}{index of first variable}
\item{jj2}{index of second variable}
\item{iprint}{flag for printing of intermediate results, including
bivariate tables for observed versus expected assuming discretized
bivariate Gaussian}
\item{prlevel}{print.level for nlm for numerical optimization}
}
\details{The intermediate bivariate tables from iprint=T are
currently not saved in the output data structure.
bprobitnllk is optimized by polychoric().
Currently prdinal2fr() and polychoric() assume that all ordinal variables have 
the same number of categories but the code could be generalized.
bprobitwPrednllk is optimized by polychoric.wPred().

polychoric0() might fail if some ordinal category has zero counts
for one or more of the ordinal variables; in this case, use polychoric()
with some preprocessing should be better.
}
\value{
d-vector for d2v(); decimal for v2d() (as way to index ordinal vectors;

bivary d-vector for d2b();

odatfr is output of ordinal2fr() (ordinal to frequencies);

a (ncat-1)xd matrix of cutpoints on the uniform(0,1) scale for unifcuts(),
where ncat is number of ordinal categories and d is number of variables;

a polychoric correlation and SE for polychoric.bivtab();

for polychoric() and polychoric0(), $polych is
a polychoric correlation matrix  based on two-stage estimate;
$iposdef is an indicator if the 2-stage correlation matrix estimate
is positive definite;

for polychoric(), also
$zcuts is the matrix of cutpoints on N(0,1) scale.

for polychoric.wPred(), $polych is
a polychoric correlation matrix  based on two-stage estimate;
$iposdef is an indicator if the 2-stage correlation matrix estimate
is positive definite
}
\seealso{
\code{\link{ordinal.bivcop}}
\code{\link{ordprobit.univar}}
}
\examples{
d2v(3,3,0) # 1 1 1 since ii=0 -> 111, ii=1 -> 112, ii=2 -> 113 etc
d2v(3,3,0,izero=TRUE) # 0 0 0 since ii=0 -> 000, ii=1 -> 001, ii=2 -> 002 etc
d2v(3,3,10) # 2 1 2 since ii=8 -> 133, next in lexicographic order are 211,212
d2v(3,2,6)  # 2 2 1  (1s and 2s)
d2v(3,2,6,izero=TRUE) # 1 1 0 (0s and 1s)
d2b(3,6)  # 1 1 0  (0s and 1s)
v2d(c(1,1,0,1),2)  # 1*2^3+ 1*2^2+ 1= 13
v2d(c(1,1,0,1),3)  # 1*3^3+ 1*3^2+ 1= 37
# examples for unifcuts and ordinal2fr
set.seed(12345)
x=rnorm(1000)
y=0.5*x+sqrt(.75)*rnorm(1000)
x=cut(x,c(-Inf,-.75,1.3,Inf))
y=cut(y,c(-Inf,-.5,1.5,Inf))
x=as.numeric(x)
y=as.numeric(y)
odat0=cbind(x,y)-1
odat1=cbind(x,y)
print(unifcuts(odat0))
print(unifcuts(odat1)) # same as above
odat0fr=ordinal2fr(odat0,3)
odat1fr=ordinal2fr(odat1,3) # same as above
# example for polychoric
data(ltmconv)
nitem=ncol(sci) # sci is 392x7
nc=4
odatfr=ordinal2fr(sci,nc) # 298x8
d=ncol(odatfr)-1 # 7
ucuts=unifcuts(sci)
zcuts=qnorm(ucuts) 
zcuts=rbind(rep(-6,d),zcuts,rep(6,d))
cat("\ncut points on N(0,1) scale\n")
print(zcuts)
polyr=polychoric(odatfr,zcuts,iprint=TRUE,prlevel=0)
polyr0=polychoric0(sci,iprint=TRUE) # without preprocessing
print(polyr0$polych-polyr$polych)
# example with bivariate table that is not square
btab=matrix(c(20,14,3,12,15,10,3,14,21,4,6,20),3,4)
polychoric.bivtab(btab,iprint=TRUE)
# example for polychoric.wPred
data(ordinalex)
xvec=c(t(ordinalex$xx))
yvec=c(t(ordinalex$yy))
ord.univar=ordprobit.univar(xvec,yvec,iprint=TRUE)
print(ord.univar)
ordtr=mord2uu(xvec,yvec,4,ord.univar$cutpts,ord.univar$beta)
polyr2=polychoric.wPred(ordtr$zzdat,iprint=TRUE)
}
\keyword{dependence}
\keyword{maximum likelihood}
\keyword{regression}
