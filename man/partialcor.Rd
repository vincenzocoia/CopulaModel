\name{partialcor}
\Rdversion{1.1}
\alias{allpcor}
\alias{partcor}
\title{
partial correlations from a correlation matrix
}
\description{
partial correlations from a correlation matrix,
(all possible) or (single)
}
\usage{
allpcor(rr)
partcor(S,given,j,k)
}
\arguments{
\item{rr}{dxd correlation matrix}
\item{S}{dxd covariance or correlation matrix}
\item{given}{the indices of the "given" or conditioning variables}
\item{j}{index of first conditioned variable}
\item{k}{index of second conditioned variable}
}
\value{
\item{pc}{single partial correlation from partcor of variables j,k given 
indices in 'given'}
\item{pcobj}{from allpcor, list with partial correlations in 2 forms: a 
3dim array named $pc3.array and a 2dim array named $pc2.array, also
$mnmx for min/max partial correlation by conditioning set size; for
pc3.array: third dimension comes from conditioning on sets in the order
1 2 12 3 13 23 123 4 14 24 124 34 134 234 1234
5 15 25 125 35 135 235 1235 45 145 245 1245 345, ...
(third dimension has length 2^d-4)
}
}
\examples{
d=5
rr=toeplitz(c(1,.5,.25,.125,.05))
pcobj=allpcor(rr)
print(pcobj$pc3.array[,,1])
partcor(rr,c(1),3,4)
partcor(rr,c(1),3,5)
print(pcobj$pc3.array[,,3])
partcor(rr,c(1,2),3,4)
partcor(rr,c(1,2),3,5)
}
\keyword{partial correlation}
