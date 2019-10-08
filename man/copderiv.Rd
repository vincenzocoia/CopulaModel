\name{copderiv}
\Rdversion{1.1}
\alias{logdfrk.deriv}
\alias{pcondfrk.deriv}
\alias{logdgum.deriv}
\alias{pcondgum.deriv}
\alias{logdbvtcop.deriv}
\alias{pcondbvtcop.deriv}
\alias{logdbb1.deriv}
\alias{pcondbb1.deriv}
\title{
Derivatives of log copula density and copula conditional cdf 
}
\description{
Derivatives of log copula density and copula conditional cdf 
}
\usage{
logdfrk.deriv(u,v,cpar)  # derivatives of c_{12}(u,v,cpar)
pcondfrk.deriv(v,u,cpar) # deriv of C_{2|1}(v+u;cpar)
# also frk can be replaced with gum or bb1 or ..
logdbvtcop.deriv(u,v,param,df=dfdefault)
pcondbvtcop.deriv(v,u,param,df=dfdefault)
}
\arguments{
\item{u}{value in interval 0,1; could be a vector}
\item{v}{value in interval 0,1; could be a vector}
\item{cpar}{copula
parameter: could be scalar or vector depending on the copula family}
\item{param}{for t copula, this is either rho in (-1,1) or
(rho,df) where df>0; in the former case, set dfdefault before using}
\item{df}{global default shape parameter for t copula}
}
\details{
These are templates. A user can add other functions like these as needed.
These could be useful for R-vine negative log-likelihoods with analytic 
derivatives for input into nlm. Result is faster with t copulas
with fixed shape parameters.
}
\value{
log pdf or \eqn{C_{2|1}} value(s) with derivatives with respect to u, v and cpar.
Output is a vector of length(3+length(cpar)) with u,v are scalars;
otherwise output is a matrix of nx [length(3+length(cpar))] if
at least one of u or v is a vector of length n.
}
\seealso{
\code{\link{pcop}}
\code{\link{pcond}}
}
\examples{
u=seq(.1,.9,.2)
v=c(.2,.2,.4,.6,.8)
dfdefault=5
logdfrk.deriv(u,v,2)
logdbvtcop.deriv(u,v,.6)
logdbvtcop.deriv(u,v,c(.6,5))
logdbb1.deriv(u,v,c(.5,2))
pcondfrk.deriv(u,v,2)
pcondbvtcop.deriv(u,v,.6)
pcondbvtcop.deriv(u,v,c(.6,5))
pcondbb1.deriv(u,v,c(.5,2))
}

\keyword{copula}
