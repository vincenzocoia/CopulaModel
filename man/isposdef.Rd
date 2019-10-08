\name{isposdef}
\Rdversion{1.1}
\alias{isposdef}
\title{
Check if Hessian matrix is positive definite
}
\description{
Check if Hessian matrix is positive definite via attempted Cholesky
decomposition
}
\usage{
isposdef(amat)
}
\arguments{
\item{amat}{dxd symmetric matrix}
}
\value{
True or False; True if matrix is positive definite
}
\examples{
a1=matrix(c(1,.5,.5,1),2,2)
a2=matrix(c(1,1.5,1.5,1),2,2)
t1=try(chol(a1))
t2=try(chol(a2))
print(isposdef(a1))
print(isposdef(a2))
}
