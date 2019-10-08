\name{factanal2nllk}
\Rdversion{1.1}
\alias{factanal2nllk}
\title{
negative log-likelihood via factor analysis
}
\description{
negative log-likelihood via factanal, assuming data are standardized
}
\usage{
factanal2nllk(rmat,mxfactor,n,iprint=F)
}
\arguments{
\item{rmat}{correlation matrix of data}
\item{mxfactor}{maximum number of factors}
\item{n}{sample size}
\item{iprint}{print flag for intermediate results}
}
\value{
vector of length mxfactor with the negative log-likelihoods for
factor models with 1,...,mxfactor factors
}
\seealso{
\code{\link{factanal.bi}}
}
\examples{
data(euro07gf);
zdat=euro07gf$zscore # 239 x 7
rr=cor(zdat)
n=nrow(zdat)
out=factanal2nllk(rr,3,n,iprint=TRUE)
}
\keyword{factor model}
