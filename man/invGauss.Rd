\name{invGauss}
\Rdversion{1.1}
\alias{invGauss}
\alias{pIG}
\alias{dIG}
\alias{qIG}
\alias{rIG}
\title{
Inverse Gaussian distribution
}
\description{
Inverse Gaussian distribution
}
\usage{
pIG(x,mu,vsi)  # cdf with mean parameter mu and second parameter vsi
dIG(x,mu,vsi)  # density
qIG(p,mu,vsi,mxiter=10,eps=1e-06,mxstep=5,iprint=F)  # inverse cdf
rIG(n,mu,vsi)  # simulation of random variables
}
\arguments{
\item{x}{positive value; could be a vector }
\item{p}{value in (0,1); could be a vector }
\item{mu}{mean parameter m of inverse Gaussian}
\item{vsi}{second parameter varsigma=lambda of inverse Gaussian}
\item{mxiter}{maximum number of iterations}
\item{eps}{tolerance for convergence}
\item{mxstep}{bound on step size for Newton-Raphson iterations}
\item{iprint}{print flag for iterations}
\item{n}{simulation sample size}
}
\value{
cdf or pdf or quantile or random sample 
}
\seealso{
\code{\link{invGaussconvfactor}}
}
\references{
Seshadri V (1993). The Inverse Gaussian Distribution. Clarendon Press.
}
\details{
Seshadri (1993): with mu>0, vsi>0 means the pdf is:
f(x;mu,vsi)= [sqrt(vsi)/sqrt(2*pi*x^3)] * exp[-(vsi/[2mu^2])*(x-mu)^2/x]
for  x>0.

Reparametrization has mu=zeta*eta as the mean and vsi=eta^2
where eta=convolution parameter;
zeta is like a scale parameter that can be set to 1 for the copula

Property of the cdf is pIG(x,mu,vsi)=pIG(x/mu,1,vsi/mu)
}
\examples{
mu=3; vsi=2
x=1:10
p=pIG(x,mu,vsi)
xx=qIG(p,mu,vsi)
print(cbind(x,p,xx))
eps=1.e-5
peps=pIG(x+eps,mu,vsi)
print(cbind((peps-p)/eps,dIG(x,mu,vsi)))
set.seed(123)
n=10000
x=rIG(n,mu,vsi)
print(summary(x))
print(var(x))
cat("theoretical values mean=", mu, " var=", mu^3/vsi,"\n")
}

\keyword{distribution}
