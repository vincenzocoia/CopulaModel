\name{btsampleStaty}
\Rdversion{1.1}
\alias{btsampleStaty}
\title{
Bootstrap from a stationary time series
}
\description{
Indexes for bootstrap from a stationary time series
}
\usage{
btsampleStaty(size,p)
}
\arguments{
\item{size}{length of time series}
\item{p}{geometric rate such as size^(-1/3)}
}
\references{
Politis DN and Romano JP (1994).
The stationary bootstrap.
Journal of the American Statistical Association,
89, 1303--1313.
}

\value{
vector of length 'size', each element is an integer between 1 and 'size'
}
\examples{
set.seed(123)
n=500
p=n^(-1/3)
for(isim in 1:5)
{ ii=btsampleStaty(n,p)
  print(length(unique(ii))) 
}
}
\keyword{resampling}
\keyword{ts}
