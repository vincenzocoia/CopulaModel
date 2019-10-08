\name{copextreme}
\Rdversion{1.1}
\alias{copextreme}
\alias{pindepcop}
\alias{pcomonocop}
\alias{pcountermono}
\alias{pcondindep}
\alias{pcondcomono}
\alias{pcondcountermono}
\title{
Boundary copula cdfs 
}
\description{
Boundary copula cdfs of independence, comonotonicity and countermonotonicity
}
\usage{
pindepcop(u,v=-1,cpar=0)
pcomonocop(u,v=-1,cpar=0)
pcountermono(u,v=-1,cpar=0)
pcondindep(v,u,cpar=0)
pcondcomono(v,u,cpar=0)
pcondcountermono(v,u,cpar=0)
}
\arguments{
\item{u}{value in interval 0,1 for bivariate; vector of values in (0,1)
for d-variate}
\item{v}{value in interval 0,1 for bivariate; -1 for d-variate}
\item{cpar}{copula parameter: could be anything but is not used;
this is for consistency when pass to functions using copula cdfs}
}
\value{
cdf value(s): pcountermono uses only u[1],u[2] if v=-1
}
\examples{
pindepcop(.1,.2)
pindepcop(seq(.1,.5,.2))
pcomonocop(.1,.2)
pcomonocop(seq(.1,.5,.2))
pcountermono(.5,.7)
pcountermono(seq(.5,.9,.2))
pcondindep(.1,.2)
pcondcomono(.1,.2)
pcondcountermono(.9,.2)
}

\keyword{copula}
