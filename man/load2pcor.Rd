\name{load2pcor}
\Rdversion{1.1}
\alias{load2pcor}
\alias{pcor2load}
\alias{grotate2}
\alias{grotate3}
\title{
Operations of loading matrix of Gaussian factor analysis
}
\description{
Operations of loading matrix of Gaussian factor analysis
}
\usage{
load2pcor(amat)   
pcor2load(rhmat)
grotate2(amat,row,iprint=F) # apply one Givens rotation
grotate3(amat,row1=1,row2=2,iprint=F) # apply three Givens rotations
}
\arguments{
\item{amat}{dxp matrix of loadings for load2pcor, dx2 for grotate2,
dx3 matrix for grotate3}
\item{rhmat}{dxp matrix for pcor2load, correlations with factor 1
in column 1, partial correlations with factor k given previous factors
in column k}
\item{row}{index of row to set to 0 in second column}
\item{row1}{index of row to set to 0 in second and third columns}
\item{row2}{index of row to set to 0 in third column}
\item{iprint}{print flag for intermediate results}
}
\value{
For load2pcor,
dxp matrix with correlations with factor 1 in the first column and
partial correlations with factor k given previous factors in the kth column.

For grotate2 and grotate3, a rotated loading matrix with same dimension
as inputted matrix.

For pcor2load, a matrix of loadings.
}
\examples{
d=7
be1=c(.7,.6,.7,.6,.7,.6,.7)
be2=c(.4,.4,.4,.4,.3,.3,.3)
cpar1.gum=gum.b2cpar(be1)
cpar2.gum=gum.b2cpar(be2)
n=300
set.seed(123)
gumdat=sim2fact(n,cpar1.gum,cpar2.gum,qcondgum,qcondgum,"gum","gum")
zdat=nscore(gumdat,iopt=TRUE)
rmat=cor(zdat)
out2=factanal(covmat=rmat,factors=2)
amat2=matrix(c(out2$loadings),d,2)
amat2=grotate2(amat2,row=1)
print(amat2) # 0 in (1,2) position
rhmat=load2pcor(amat2)
amat=pcor2load(rhmat)
out3=factanal(covmat=rmat,factors=3)
amat3=matrix(c(out3$loadings),d,3)
amat3=grotate3(amat3,row1=1,row2=2)
print(amat3) # 0 in (1,2), (1,3) and (2,3) positions
load2pcor(amat3)
}
\keyword{factor model}
\keyword{partial correlation}
