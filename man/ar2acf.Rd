\name{ar2acf}
\Rdversion{1.1}
\alias{ar2acf}
\title{
Autocorrelation of Gaussian AR2 given lag1 and lag2 correlations
}
\description{
Autocorrelation of Gaussian AR2 given lag1 and lag2 correlations
}
\usage{
ar2acf(rho1,rho2,d)
}
\arguments{
\item{rho1}{correlation for lag 1}
\item{rho2}{correlation for lag 2, so correlation matrix with
rho1,rho1,rho2 should be positive definite}
\item{d}{dimension of desired Toeplitz matrix}
}
\value{
acf function up to lag d-1
}
\examples{
acf10=ar2acf(.6,.5,10)
toeplitz(acf10[1:6])
}
\keyword{ts}
