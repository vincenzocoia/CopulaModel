\name{invGaussconvfactor}
\Rdversion{1.1}
\alias{invGaussconvfactor}
\alias{rIGconv}
\alias{pmIGfact}
\alias{pmIGfact.gl}
\alias{dmIGfact}
\alias{dmIGfact.gl}
\alias{dbIGfact}
\alias{pmIGfcop.gl}
\alias{dmIGfcop.gl}
\title{
Inverse Gaussian convolution factor model
}
\description{
Inverse Gaussian convolution factor model
}
\usage{
rIGconv(n,th0,thvec,ze=1)
pmIGfact(xvec,th0,thvec,zero=0)
pmIGfact.gl(xvec,th0,thvec,gl)
dmIGfact(xvec,th0,thvec,zero=0)
dmIGfact.gl(xvec,th0,thvec,gl)
dbIGfact(x1,x2,th0,th1,th2,zero=0)  
pmIGfcop.gl(uvec,param,gl)
dmIGfcop.gl(uvec,param,gl)
}
\arguments{
\item{n}{sample size for simulation}
\item{th0}{scalar for shape parameter of the shared/common component}
\item{thvec}{vector of shape parameters of individual components, length d}
\item{param}{parameter vector with length d+1 with th0,thvec}
\item{xvec}{vector of length d with positive values}
\item{uvec}{vector of length d with values in (0,1)}
\item{gl}{Gauss-Legendre object with components $nodes and $weights}
\item{ze}{non-convolution parameter zeta, can be set to 1 for copula}
\item{zero}{tolerance for numerical integration, set as 0.0001 if there
are numerial problems}
\item{x1}{positive value for first variable (bivariate case)}
\item{x2}{positive value for second variable (bivariate case)}
\item{th1}{scalar for shape parameter of first variable (bivariate case)}
\item{th2}{scalar for shape parameter of second variable (bivariate case)}
}
\value{
random sample (nxd matrix) for rIGconv

cdf or pdf for remaining functions
}
\seealso{
\code{\link{invGauss}}
\code{\link{gammaconvfactor}}
}
\examples{
n=1000
th0=2
thvec=c(.3,.3)
set.seed(123)
xdat=rIGconv(n,th0,thvec)
cor(xdat)
#plot(xdat)
gl=gausslegendre(25)
pmIGfact(c(1,1.1),th0,c(.4,.4),zero=0)
pmIGfact.gl(c(1,1.1),th0,c(.4,.4),gl)
# check that density is finite on diagonal
dbIGfact(1,1.1,th0,th1=.4,th2=.4,zero=0)
dmIGfact(c(1,1.1),th0,c(.4,.4),zero=0)
dmIGfact.gl(c(1,1.1),th0,c(.4,.4),gl)
dbIGfact(1,1.0001,th0,th1=.4,th2=.7,zero=0)
# copula 
pmIGfcop.gl(c(.5,.6),c(2,.4,.4),gl)
dmIGfcop.gl(c(.5,.6),c(2,.4,.4),gl)
}
\keyword{factor model}
