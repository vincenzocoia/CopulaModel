\name{rectmult}
\Rdversion{1.1}
\alias{rectmult}
\alias{rectmfrk}
\alias{rectmgum}
\alias{rectmgal}
\alias{rectemvn}
\title{
multivariate rectangle probabilities for copula models
}
\description{
multivariate rectangle probabilities for copula models with discrete
responses
}
\usage{
rectmult(u1vec,u2vec,cpar,pmcop)
rectmfrk(u1vec,u2vec,cpar)  # multivariate exchangeable Frank
rectmgum(u1vec,u2vec,cpar)  # multivariate exchangeable Gumbel
rectmgal(u1vec,u2vec,cpar)  # multivariate exchangeable Galambos
rectemvn(u1vec,u2vec,cpar)  # multivariate exchangeable normal
}
\arguments{
\item{u1vec}{vector of dimension d, lower corner of rectangle}
\item{u2vec}{vector of dimension d, upper corner of rectangle}
\item{cpar}{parameter of the d-variate copula}
\item{pmcop}{function for the cdf of the d-variate copula}
}
\value{rectangle probability
}
\seealso{
\code{\link{copmultiv}}
}
\examples{
rectmgum(c(.1,.2,.1),c(.8,.9,.7),2)
rectmgum(c(.1,.2,.1),c(.8,.9,.7),4)
rectmult(c(.1,.2,.1),c(.8,.9,.7),4,pmcop=pmgum)
rectmfrk(c(.1,.2,.1),c(.8,.9,.7),2)
rectmfrk(c(.1,.2,.1),c(.8,.9,.7),4)
rectmgal(c(.1,.2,.1),c(.8,.9,.7),3)
rectemvn(c(.1,.2,.1),c(.8,.9,.7),0.6)
}
\keyword{copula}
\keyword{distribution}
