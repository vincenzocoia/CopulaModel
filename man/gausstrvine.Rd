\name{gausstrvine}
\Rdversion{1.1}
\alias{gausstrvine}
\alias{gausstrvine.nonuniq}
\title{
Best Gaussian truncated vines
}
\description{
Best Gaussian truncated d-dimensional vines up to d-2 trees
}
\usage{
gausstrvine(rmat,iprint=F)
gausstrvine.nonuniq(rmat,jtrunc=3,eps=1.e-7,iprint=F) # check non-uniqueness
}
\arguments{
\item{rmat}{dxd correlation matrix, 4<=d<=8}
\item{iprint}{print flag for intermediate steps (in f90 code)}
\item{jtrunc}{truncation level to check on degree of non-uniqueness}
\item{eps}{tolerance to check on degree of non-uniqueness, default 1e-7}
}
\value{
\item{bnum}{d-2 dimensional vector with indices of best vine arrays;
can get vine arrays with something like vnum2array(d, bnum[ell])}
\item{logdetmx}{d-1 dimensional vector with max log determinants
for truncated vines of order 1,...,d-2; the last entry in position
d-1 is the log determinant of rmat}
\item{permmat}{dx(d-2) matrix, with permutation leading to a
best ell-truncated vine in column ell}
\item{pcarr}{dxdx(d-2) partial correlation array with matrix of
partial correlations in pcarr[,,ell] for the best ell-truncated vine}
}
\seealso{
\code{\link{gausstrvineMST}}
}
\details{
Note that even if the optimal ell-truncated vine is unique, the
vine array leading to it is not unique.
The output of gausstrvine.nonuniq() is not saved into R variables,
so it should be looked at for other truncated vines that lead to the same
determinant as the optimal.
}
\examples{
\dontrun{
rmat=matrix(c(
1.00000,0.69965,0.70477,0.66536,0.65967,
0.69965,1.00000,0.65499,0.61713,0.61202,
0.70477,0.65499,1.00000,0.62967,0.62798,
0.66536,0.61713,0.62967,1.00000,0.57398,
0.65967,0.61202,0.62798,0.57398,1.00000), 5,5)
out=gausstrvine(rmat,iprint=FALSE)
print(out)
outnonuniq=gausstrvine.nonuniq(rmat,jtrunc=3,eps=1.e-7,iprint=TRUE)
# some checks
d=nrow(rmat)
for(ell in 1:(d-2))
{ A=vnum2array(d,out$bnum[ell])
  cat("truncation level ", ell,"\n")
  print(A)
  cat("check on log determinant\n")
  pcmat=out$pcarr[,,ell]
  logdet=sum(log(1-pcmat[1:ell,]^2))
  print(logdet)
}
print(determinant(rmat,log=TRUE)$modulus) }
}
\keyword{partial correlation}
\keyword{vine}
