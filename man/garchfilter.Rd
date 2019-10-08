\name{garchfilter}
\alias{garchfilter}
\alias{gfiltersubset}
\title{GARCH filter applied separately to log returns}
\description{
GARCH filter applied separately to a matrix of log returns, 
the d log returns should have the same dates, a subset of rows can be chosen 
}
\usage{
gfiltersubset(lgret,ar,m1,m2,iprint=F)
}
\arguments{
\item{lgret}{Nxd matrix of log returns of d financial assets}
\item{ar}{TRUE for univariate GARCH(1,1)-AR(1) for model, 
FALSE for GARCH(1,1) model}
\item{m1}{first row of lgret to use for the subset}
\item{m2}{last row of lgret to use for the subset}
\item{iprint}{print flag for GARCH estimates, F by default}
}
\value{
\item{filter}{GARCH filtered data (nxd, where n=m2-m1+1)}
\item{uscore}{empirical uniform scores (nxd)}
\item{zscore}{empirical normal scores (nxd)}
\item{uscmodel}{model-based uniform scores (nxd)}
\item{zscmodel}{model-based normal scores (nxd)}
\item{sigmat}{matrix of estimated volatilities (nxd)}
\item{coef}{matrix of GARCH parameters (6xd or 5xd for ar=T or ar=F
respectively), the parameters are mu, (ar1), omega, alpha1, beta1, shape.}
}
\details{
This function requires library fGarch and its dependent libraries.
}
\examples{
\dontrun{
data(euro07lr) # euro07lr data set is 239x7
library(fGarch)
out=gfiltersubset(euro07lr,ar=FALSE,1,200,iprint=TRUE) }
}
\keyword{ts}
