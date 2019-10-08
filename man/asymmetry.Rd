\name{asymmetry}
\Rdversion{1.1}
\alias{asymmetry}
\alias{skewrefl}
\alias{skewperm}
\alias{qskewrefl}
\alias{qskewperm}
\title{
Asymmetry / skewness measures for bivariate copulas
}
\description{
Reflection and permutation asymmetry / skewness measures for bivariate copulas
}
\usage{
skewrefl(uu)  # third moment (u1+u2-1)^3
skewperm(uu)  # third moment  (u1-u2)^3
qskewperm(uu,p=.05,nrep=100)  # quantile-based measure for u1-u2
qskewrefl(uu,p=.05,nrep=100)  # quantile-based measure for u1+u2-1
}
\arguments{
\item{uu}{nx2 data matrix of uniform scores}
\item{p}{value between 0 and 0.5, leading to quantiles at p and 1-p}
\item{nrep}{number of replications for bootstrap SE}
}
\value{skewness/asymmetry measure, with a crude SE
}
\references{
Rosco J-F and Joe H (2014). Measures of tail asymmetry for bivariate
copulas. Statistical Papers, 54, 709-726.
}
\examples{
n=100
set.seed(123)
urefl=rbreflasym(n,.25)  # extreme reflection asymmetry
uperm=rbpermasym(n,.25)  # extreme permutation asymmetry
n=1000
set.seed(123)
uagum=rbasymgum1(n,c(2,.4))  # bivariate asymmetric Gumbel
umo=rbMO1(n,p=.4) # bivariate 1-parameter Marshall-Olkin
\dontrun{
par(mfrow=c(2,2))
plot(urefl)
plot(uperm)
plot(uagum)
plot(umo) }
sk1=skewrefl(urefl)
sk2=skewperm(uperm)
skagum=skewperm(uagum)
skmo=skewperm(umo)
print(rbind(sk1,sk2,skagum,skmo))
qsk1=qskewrefl(urefl)
qsk2=qskewperm(uperm)
qskagum=qskewperm(uagum)
qskmo=qskewperm(umo)
print(rbind(qsk1,qsk2,qskagum,qskmo)) }
\keyword{nonparametric}
