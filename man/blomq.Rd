\name{blomq}
\Rdversion{1.1}
\alias{blomq}
\title{
Blomqvist's beta for bivariate data set
}
\description{
Blomqvist's beta for bivariate data set
}
\usage{
blomq(bivdat,iunif=F)
}
\arguments{
\item{bivdat}{data matrix with 2 columns }
\item{iunif}{flag for whether bivdat has already been converted to
uniform(0,1) scores, default=F in which case conversion is done}
}
\value{
empirical Blomqvist's beta
}
\examples{
set.seed(1234)
n=300; cpar=3
udat=rgum(n,cpar)
blomq(udat)
4*pgum(.5,.5,cpar)-1  # theoretical value
}
\keyword{dependence}
