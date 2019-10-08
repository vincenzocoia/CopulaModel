\name{cormat}
\Rdversion{1.1}
\alias{corDis}
\alias{corvec2mat}
\alias{cormat2vec}
\title{
Functions on correlation matrices
}
\description{
Functions on correlation matrices
}
\usage{
corvec2mat(rvec) # convert correlation vector to matrix
cormat2vec(rmat) # extract vector from correlation matrix
corDis(Rmod,Robs,n=0,npar=0) # discrepancy between model-based and observed correlation matrices
}
\arguments{
\item{rvec}{vector of length d*(d-1)/2 with order 
r[1,2],r[1,3],r[2,3],r[1,4], ..., r[d-1,d]}
\item{rmat}{dxd correlation matrix}
\item{Rmod}{model-based correlation matrix}
\item{Robs}{observed correlation matrix}
\item{n}{sample size for Robs}
\item{npar}{parameter vector size leading to Rmod}
}
\value{
dxd correlation matrix for corvec2mat

vector of length d*(d-1)/2 for cormat2vec

log(det(Rmod))-log(det(Robs))+sum(diag(solve(Rmod,Robs)))-nrow(Robs)
for corDis assuming a Gaussian dependence model
}
\examples{
rvec=c(.3,.4,.5,.4,.6,.7)
Rmod=corvec2mat(rvec)
print(Rmod); print(chol(Rmod))
print(cormat2vec(Rmod))
robsvec=c(.32,.38,.53,.41,.61,.67)
Robs=corvec2mat(robsvec)
print(corDis(Rmod,Robs))
print(corDis(Rmod,Robs,n=400,npar=3))
}
\keyword{dependence}
