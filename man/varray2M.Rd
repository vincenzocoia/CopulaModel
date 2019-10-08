\name{varray2M}
\Rdversion{1.1}
\alias{varray2M}
\title{
Vine array to maximum matrix
}
\description{
Vine array to maximum matrix
}
\usage{
varray2M(A,iprint=F,str="")
}
\arguments{
\item{A}{dxd vine array with 1:d on diagonal; only upper triangle is used}
\item{iprint}{print flag for some intermediate steps}
\item{str}{string to describe the vine if iprint=T }
}
\value{
\item{mxarray}{array with M[k,j]= max A[k,1],..,A[k,j]}
\item{icomp}{array with indicators for use in log-likelihood and simulation}
}
\details{In the VineCopula R package and in 
"Dissmann J, Brechmann EC, Czado C and Kurowicka D (2013).
Computational Statistics and Data Analysis, 59, 52-69",
the vine array is abbreviation as RVM and it is A[d:1,d:1]
where d=ncol(A). This conversion is needed if this package and
VineCopula are used together.
}
\seealso{
\code{\link{vinearray}}
}
\examples{
C= matrix(c(1,0,0,0,0, 1,2,0,0,0, 1,2,3,0,0, 1,2,3,4,0, 1,2,3,4,5), 5,5)
D= matrix(c(1,0,0,0,0, 1,2,0,0,0, 2,1,3,0,0, 3,2,1,4,0, 4,3,2,1,5), 5,5)
B0=matrix(c(1,0,0,0,0, 1,2,0,0,0, 1,2,3,0,0, 1,2,3,4,0, 1,3,2,4,5), 5,5)
B1=matrix(c(1,0,0,0,0, 1,2,0,0,0, 1,2,3,0,0, 1,2,3,4,0, 2,1,3,4,5), 5,5)
B2=matrix(c(1,0,0,0,0, 1,2,0,0,0, 2,1,3,0,0, 1,2,3,4,0, 1,2,3,4,5), 5,5)
B3=matrix(c(1,0,0,0,0, 1,2,0,0,0, 1,2,3,0,0, 1,3,2,4,0, 2,1,3,4,5), 5,5)
D6= matrix(c(1,0,0,0,0,0, 1,2,0,0,0,0, 2,1,3,0,0,0, 3,2,1,4,0,0, 4,3,2,1,5,0,
   5,4,3,2,1,6),6,6)
varray2M(C ,iprint=TRUE,"C ")
varray2M(D ,iprint=TRUE,"D ")
varray2M(B0,iprint=TRUE,"B0")
varray2M(B1,iprint=TRUE,"B1")
varray2M(B2,iprint=TRUE,"B2")
varray2M(B3,iprint=TRUE,"B3")
varray2M(D6,iprint=TRUE,"D6")
}
\keyword{vine}
