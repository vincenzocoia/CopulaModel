\name{bivpmf}
\Rdversion{1.1}
\alias{bivpmf2cdf}
\alias{corbivpmf}
\title{
bivariate probability mass function: cdf and correlation  
}
\description{
bivariate probability mass function: cdf and correlation    
}
\usage{
bivpmf2cdf(bpmf)
corbivpmf(bpmf)
}
\arguments{
\item{bpmf}{bivariate probability mass function in a matrix
(non-negative and sums to 1)}
}
\value{
matrix with bivariate cdf for bivpmf2cdf()

means, variances, covariance, correlation for corbivpmf(),
assuming categories are 0:(nrow(pbmf)-1) and
 0:(ncol(pbmf)-1)
}
\examples{
bpmf=matrix(c(.1,.2,.3,.2,.1,.1),2,3)
bivpmf2cdf(bpmf)
corbivpmf(bpmf)
}
\keyword{distribution}
