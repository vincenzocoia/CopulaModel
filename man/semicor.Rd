\name{semicor}
\Rdversion{1.1}
\alias{semicor}
\alias{semicortable}
\title{
Semi-correlations for bivariate or multivariate data set
}
\description{
Semi-correlations for bivariate or multivariate data set
}
\usage{
semicor(bivdat,inscore=T)
semicortable(mdat,inscore=F)
}
\arguments{
\item{bivdat}{data matrix with 2 columns }
\item{mdat}{data matrix with d columns }
\item{inscore}{T if bivdat or mdat has already been converted to normal scores}
}
\value{
For semicor,
$ncorr: correlation of normal scores;
$lcorr: lower semi-correlation of normal scores;
$ucorr: upper semi-correlation of normal scores.

For semicortable, d*(d-1)/2 by 6 table with 
columns j1,j2,ncorr,lcorr,ucorr,bvnsemic
}

\seealso{
\code{\link{bvnsemic}}
\code{\link{rhoNsemic}}
\code{\link{tailweightedDepmeas}}
}
\examples{
set.seed(1234)
n=300
x=matrix(rnorm(2*n),n,2)
rho=0.6
x[,2]=rho*x[,1]+sqrt(1-rho^2)*x[,2]
semicor(x,inscore=TRUE)
semicor(x,inscore=FALSE)
bvnsemic(rho)  # theoretical value
# log return data
data(euro0306)
semicortable(euro0306[,2:6],inscore=FALSE)
}
\keyword{dependence}
