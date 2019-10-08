\name{pcond}
\Rdversion{1.1}
\alias{pcond}
\alias{qcond}
\alias{pcondbvncop}
\alias{pcondbvtcop}
\alias{pcondt}
\alias{pcondpla}
\alias{pcondfrk}
\alias{pcondmtcj}
\alias{pcondmtcjr}
\alias{pcondjoe}
\alias{pcondgum}
\alias{pcondgumr}
\alias{pcondgal}
\alias{pcondhr}
\alias{pcondfgm}
\alias{pcondbMO21}
\alias{pcondbb1}
\alias{pcondbb1r}
\alias{pcondbb2}
\alias{pcondbb3}
\alias{pcondbb4}
\alias{pcondbb5}
\alias{pcondbb6}
\alias{pcondbb7}
\alias{pcondbb7r}
\alias{pcondbb8}
\alias{pcondbb9}
\alias{pcondbb10}
\alias{pcondtev}
\alias{qcondbvncop}
\alias{qcondbvtcop}
\alias{qcondt}
\alias{qcondpla}
\alias{qcondfrk}
\alias{qcondmtcj}
\alias{qcondmtcjr}
\alias{qcondjoe}
\alias{qcondgum}
\alias{qcondgumr}
\alias{qcondgal}
\alias{qcondhr}
\alias{qcondfgm}
\alias{qcondbb1}
\alias{qcondbb1r}
\alias{qcondbb2}
\alias{qcondbb3}
\alias{qcondbb4}
\alias{qcondbb5}
\alias{qcondbb6}
\alias{qcondbb7}
\alias{qcondbb7r}
\alias{qcondbb8}
\alias{qcondbb9}
\alias{qcondbb10}
\title{
Bivariate copula conditional cdfs and quantile functions
}
\description{
Bivariate copula conditional cdfs and quantile functions
}
\usage{
pcond(v,u,cpar)
qcond(p,u,cpar)
}
\arguments{
\item{v}{conditioned value in interval 0,1; could be a vector }
\item{u}{conditioning value in interval 0,1; could be a vector }
\item{p}{quantile in interval 0,1; could be a vector }
\item{cpar}{copula parameter: could be scalar or vector depending on the 
copula family, could be a matrix with m columns if copula family has
m parameters. }
}
\details{
Choices appending 'pcond' and 'qcond'
are pla, frk, mtcj, mtcjr (reflected mtcj), joe, gum,
gumr, gal, hr, fgm, bMO21,
bb1, bb2, bb3, bb4, bb4r, bb5, bb6, bb7, bb7r, bb8, bb9,
bb10, tev, ipsA, ipsAr, imitlefA, imitlefAr, etc.
Use pcondbvncop and pcondbvtcop (or pcondt) for the conditional 
cdfs of the bivariate normal and t copulas.  

See help page for pcop for the abbreviations of the copula names.
}
\value{
conditional cdf value(s) or inverse conditional cdf value(s)
}
\references{
Joe H (1997). Multivariate Models and Dependence Concepts. Chapman & Hall.
}
\seealso{
\code{\link{pcond}}
\code{\link{cparbound}}
}
\examples{
u=seq(.1,.9,.2)
v=u
pcondpla(u,v,2)
pcondfrk(u,v,2)
pcondmtcj(u,v,2)
pcondbb1(u,v,c(.5,1.2))
qcondpla(u,v,2)
qcondfrk(u,v,2)
qcondmtcj(u,v,2)
qcondbb1(u,v,c(.5,1.2))
qcondjoe(u,v,2)
qcondgum(u,v,2)
qcondgal(u,v,2)
qcondhr(u,v,2)
}

\keyword{copula}
