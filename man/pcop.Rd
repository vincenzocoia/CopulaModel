\name{pcop}
\Rdversion{1.1}
\alias{pcop}
\alias{dcop}
\alias{logdcop}
\alias{pbvncop}
\alias{pbvtcop}
\alias{ppla}
\alias{pfrk}
\alias{pmtcj}
\alias{pmtcjr}
\alias{pjoe}
\alias{pgum}
\alias{pgumr}
\alias{pgal}
\alias{phr}
\alias{pfgm}
\alias{pbMO}
\alias{pbasymgum1}
\alias{pbb1}
\alias{pbb1r}
\alias{pbb2}
\alias{pbb3}
\alias{pbb4}
\alias{pbb4r}
\alias{pbb5}
\alias{pbb6}
\alias{pbb7}
\alias{pbb7r}
\alias{pbb8}
\alias{pbb9}
\alias{pbb10}
\alias{ptev}
\alias{dbvncop}
\alias{dbvtcop}
\alias{dpla}
\alias{dfrk}
\alias{dmtcj}
\alias{dmtcjr}
\alias{djoe}
\alias{dgum}
\alias{dgumr}
\alias{dgal}
\alias{dhr}
\alias{dfgm}
\alias{dbb1}
\alias{dbb1r}
\alias{dbb2}
\alias{dbb3}
\alias{dbb4}
\alias{dbb4r}
\alias{dbb5}
\alias{dbb6}
\alias{dbb7}
\alias{dbb7r}
\alias{dbb8}
\alias{dbb9}
\alias{dbb10}
\alias{dtev}
\alias{logdbvncop}
\alias{logdbvtcop}
\alias{logdtcop}
\alias{logdpla}
\alias{logdfrk}
\alias{logdmtcj}
\alias{logdmtcjr}
\alias{logdjoe}
\alias{logdgum}
\alias{logdgal}
\alias{logdhr}
\alias{logdfgm}
\alias{logdbb1}
\alias{logdbb1r}
\alias{logdbb7}
\title{
Bivariate copula cdfs and densities
}
\description{
Bivariate copula cdfs and densities, for parametric families
Log of copula densities.
}
\usage{
pcop(u,v,cpar)
dcop(u,v,cpar)
logdcop(u,v,cpar)
}
\arguments{
\item{u}{value in interval 0,1; could be a vector}
\item{v}{value in interval 0,1; could be a vector}
\item{cpar}{copula parameter: could be scalar or vector depending on the 
copula family, could be a matrix with m columns if copula family has
m parameters.}
}
\details{
Choices are 'cop' in pcop and dcop are bvncop, bvtcop,
pla, frk, mtcj, mtcjr (reflected mtcj), joe, gum, gumr, gal, hr, fgm, bMO,
tev, bb1, bb1r, bb2, bb3, bb4, bb4r, bb5, bb6, bb7, bb7r, bb8, bb9,
bb10, ipsA, ipsAr, imitlefA, imitlefAr, etc.
Use dbvncop and dbvtcop for the bivariate normal and t copula densities.
Note that pbvtcop assumes that the degree of freedom parameter is a
positive integer.

The bounds for the copula parameter(s) are in the source R files,
or can get obtained from the function cparbound().
The copula names are abbreviations for:

bvn = bivariate normal or Gaussian

bvt = bivariate t

pla = Plackett

frk = Frank

mtcj = Mardia-Takahasi-Clayton-Cook-Johnson

joe = Joe/B5

gum = Gumbel

gal = Galambos

hr = Huesler-Reiss

fgm = Farlie-Gumbel-Morgenstern

bMO = bivariate Marshall-Olkin (cdf only, it has not absolutely continuous)

basymgum1 = bivariate asymmetric Gumber with one skew parameter (cdf only)

tev = t-EV = extreme value limit of bivariate t

bb1 = BB1 etc

ipsA = Archimedean based on integrated positive stable LT

imitlefA = Archimedean based on integrated Mittag-Leffler LT

Choices are 'cop' in logdcop are pla, frk, mtcj, mtcjr, joe, gum, gal, hr, fgm,
bvncop, bvtcop, bb1, bb1r, bb7, ipsA, ipsAr.
These are included if it is more efficient to code logdcop directly,
Otherwise create your own function from log of the appropriate dcop
function. Another possibility is to write dcop functions with log=T option. 
}
\value{
cdf or pdf or log pdf value(s)
}
\references{
Joe H (1997). Multivariate Models and Dependence Concepts. Chapman & Hall.
}
\seealso{
\code{\link{pcond}}
\code{\link{cparbound}}
}
\examples{
u=seq(.1,.9,.1)
v=u
ppla(u,v,2)
pfrk(u,v,2)
pmtcj(u,v,2)
pbb1(u,v,c(.5,1.2))
dpla(u,v,2)
logdpla(u,v,2)
log(dpla(u,v,2))
}
\keyword{copula}
\keyword{distribution}
