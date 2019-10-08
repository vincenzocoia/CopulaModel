\name{gammaconvfactor}
\Rdversion{1.1}
\alias{gammaconvfactor}
\alias{rgammaconv}
\alias{pmgamfact}
\alias{dmgamfact}
\alias{dmgamfact.gl}
\alias{pbgamfact}
\alias{dbgamfact}
\alias{dbgamfact.gl}
\alias{pbgamfcop}
\alias{dbgamfcop}
\alias{dbgamfcop.gl}
\alias{pcondbgamfcop21}
\alias{pcondbgamfcop12}
\alias{pmgamfcop}
\alias{dmgamfcop}
\alias{dmgamfcop.gl}
\title{
Gamma convolution factor model
}
\description{
Gamma convolution factor model and its copula
}
\usage{
rgammaconv(n,th0,thvec)
pmgamfact(xvec,th0,thvec,zero=0)
dmgamfact(xvec,th0,thvec,zero=0)
dmgamfact.gl(xvec,th0,thvec,gl)
pbgamfact(x1,x2,th0,th1,th2,zero=0) 
dbgamfact(x1,x2,th0,th1,th2,zero=0)
dbgamfact.gl(x1,x2,th0,th1,th2,gl)
pbgamfcop(u,v,param)
dbgamfcop(u,v,param,zero=0)
dbgamfcop.gl(u,v,param,gl=gldefault)
pcondbgamfcop21(v,u,param,zero=0)
pcondbgamfcop12(u,v,param,zero=0)
pmgamfcop(uvec,param)
dmgamfcop(uvec,param,zero=0)
dmgamfcop.gl(uvec,param,gl)
}
\arguments{
\item{n}{sample size for simulation}
\item{th0}{scalar for shape parameter of the shared/common component}
\item{thvec}{vector of shape parameters of individual components, length d}
\item{xvec}{vector of length d with positive values}
\item{uvec}{vector of length d with values in (0,1)}
\item{gl}{Gauss-Legendre object with components $nodes and $weights}
\item{x1}{positive value for first variable (bivariate case)}
\item{x2}{positive value for second variable (bivariate case)}
\item{th1}{scalar for shape parameter of first variable (bivariate case)}
\item{th2}{scalar for shape parameter of second variable (bivariate case)}
\item{zero}{tolerance for numerical integration, set as 0.0001 if there
are problems}
\item{u}{value in interval 0,1; could be a vector }
\item{v}{value in interval 0,1; could be a vector }
\item{param}{parameter vector with length d+1 with th0,thvec}
}
\value{
random sample (nxd matrix) for rgammaconv

cdf or conditional cdf or pdf for remaining functions
}
\seealso{
\code{\link{invGaussconvfactor}}
}
\examples{
n=1000
th0=2
thvec=c(.3,.3)
set.seed(123)
xdat=rgammaconv(n,th0,thvec)
cor(xdat)
#plot(xdat)
gl=gausslegendre(35)
pmgamfact(c(1,1.1),th0,c(.4,.4),zero=0)
dbgamfact(1,1.1,th0,th1=.4,th2=.4,zero=0)
dbgamfact.gl(1,1.1,th0,th1=.4,th2=.4,gl) # could be inaccurate for th1,th2<1
dbgamfact(1,1.1,th0,th1=1.2,th2=1.4,zero=0)
dbgamfact.gl(1,1.1,th0,th1=1.2,th2=1.4,gl) 
dmgamfact(c(1,1.1),th0,c(1.2,1.4),zero=0)
dmgamfact.gl(c(1,1.1),th0,c(1.2,1.4),gl)
# density can be finite on diagonal
try(dbgamfact(1,1.0001,0.1,th1=.4,th2=.4,zero=0))
try(dbgamfact(1,1,0.1,th1=.4,th2=.4,zero=0))
# copula 
pmgamfcop(c(.5,.6),c(2,1.2,1.4))
dmgamfcop(c(.5,.6),c(2,1.2,1.4),zero=0)
dmgamfcop.gl(c(.5,.6),c(2,1.2,1.4),gl)
pcondbgamfcop21(.6,.5,c(2,1.2,1.4),zero=0)
pcondbgamfcop12(.5,.6,c(2,1.2,1.4),zero=0)
pcondbgamfcop21(.5,.6,c(2,1.4,1.2),zero=0)
}

\keyword{copula}
\keyword{factor model}
