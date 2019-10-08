\name{rLTstochrep}
\Rdversion{1.1}
\alias{rLTstochrep}
\alias{rmcop}
\alias{rpostable}
\alias{rsibuya}
\alias{rlogseries}
\alias{rmitlef}
\alias{rgammaSgamma}
\alias{rpostableSgamma}
\alias{rsibuyaSpostable}
\alias{rsibuyaSgamma}
\alias{rmfrk}
\alias{rmfrk0}
\alias{rmmtcj}
\alias{rmjoe}
\alias{rmgum}
\alias{rmbb1}
\alias{rmbb2}
\alias{rmbb3}
\alias{rmbb6}
\alias{rmbb7}
\alias{rmbb10}
\title{
Simulation of random variables from LTs used in Archimedean copulas.
Simulation from multivariate Archimedean copulas.
}
\description{
Simulation of random variables from LTs used in Archimedean copulas.
Simulation from multivariate Archimedean copulas.
}
\usage{
rpostable(n,alp) # 0<alp<1
rsibuya(n,alp) # 0<alp<1
rlogseries(n,cpar) # cpar as for Frank leads to better parametrization
rmitlef(n,param)
rgammaSgamma(n,param)
rpostableSgamma(n,param)
rsibuyaSpostable(n,param)
rsibuyaSgamma(n,param)
rmfrk0(n,d,cpar) # R version
rmfrk(n,d,cpar,icheck=F) # link to C
rmcop(n,d,cpar)  # choices for 'cop' for mtcj,joe,gum,bb1,bb2,bb3,bb6,bb7,bb10
}
\arguments{
\item{n}{sample size }
\item{d}{dimension for multivariate Archimedean}
\item{alp}{parameter of Laplace transform (LT)}
\item{param}{(vector) parameter of Laplace transform (LT)}
\item{cpar}{copula parameter: 
could be scalar or vector depending on the copula family}
\item{icheck}{flag to print out means and correlation as checks}
}
\seealso{
\code{\link{rcop}}
\code{\link{rbivcop2param}}
}
\details{
The LT families matching the Archimedean copula families are:

logseries for Frank;

gamma for MTCJ=Mardia-Takahasi-Cook-Johnson;

Sibuya for Joe;

positive stable for Gumbel;

Mittag-Leffler (or gamma stopped positive stable) for multivariate version of BB1;

gammaSgamma (gamma stopped gamma) for multivariate version of BB2;

postableSgamma (positive stable stopped gamma) for multivariate version of BB3;

sibuyaSpostable (Sibuya stopped positive stable) for multivariate version of BB6;

sibuyaSgamma (Sibuya stopped gamma) for multivariate version of BB7;

shifted negative binomial (see code) for multivariate version of BB10.
}
\value{
vector for rpostable to rsibuyaSgamma,
nxd matrix for rmfrk to rmbb10.
}
\references{
Joe H (2014). Dependence Modeling with Copulas. Chapman&Hall/CRC.
See Appendix for the names of some of the LTs and the source of the algorithms.
}
\examples{
cpar=c(1,2); n=1000
r=rmitlef(n,cpar)
print(summary(r))
uu=rmbb1(n,d=3,cpar)
print(summary(uu))
print(taucor(uu[,1],uu[,2]))
print(taucor(uu[,1],uu[,3]))
print(taucor(uu[,2],uu[,3]))
tau=bb1.cpar2tau(cpar)
cat("theor.tau=",tau,"\n")
}
\keyword{copula}
\keyword{distribution}
\keyword{simulation}
