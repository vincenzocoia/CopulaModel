\name{pbnorm}
\alias{pbnorm}
\alias{pbvt}
\title{
Bivariate normal and Student cdfs with vectorized inputs  
}
\description{
Bivariate normal and Student cdfs with vectorized inputs
}
\usage{
pbnorm(z1,z2,rho,icheck=F)
pbvt(z1,z2,param,icheck=F)
}
\arguments{
\item{z1}{scalar or vector of reals }
\item{z2}{scalar or vector of reals }
\item{rho}{scalar or vector parameter in (-1,1);
vectors cannot have different lengths if larger than 1, each of z1,z2,rho
either has length 1 or a constant n greater than 1 }
\item{param}{vector of length 2, or matrix with 2 columns;
vectors and number of rows of matrix cannot be different if larger than 1;
for param, first column is rho, second column is df.}
\item{icheck}{T if checks are made for proper inputs, default of F}
}
\value{
cdf value(s)
}
\details{Donnelly's code can be inaccurate in the tail when the tail
probability is 2.e-9 or less (it sometimes returns 0).
In the case the exchmvn code is used with dimension 2.
Alternatively a user can use vectorized function pbivnorm() in the library 
pbivnorm, and write a function pbvncop based on it.
}
\references{
Donnelly TG (1973). Algorithm 462: bivariate normal distribution, 
Communications of the Association for Computing Machinery, 16, 638;

pbvt is taken from the source of the mvtnorm R package, with a different
interface to R; the degree of freedom parameter for the multivariate t cdf
is integer value in that source. One could use linear interpolation
for non-integer-valued degree of freedom parameter.
}
\seealso{
\code{\link{exchmvn}}
}
\examples{
cat("\n pbnorm rho changing\n")
z1=.3; z2=.4; rho=seq(-1,1,.1)
out1=pbnorm(z1,z2,rho)
print(cbind(rho,out1))

cat("\n pbnorm matrix inputs for z1, z2\n")
rho=.4
z1=c(-.5,.5,10.)
z2=c(-.4,.6,10.)
z1=matrix(z1,3,3)
z2=matrix(z2,3,3,byrow=TRUE)
out3=pbnorm(z1,z2,rho)
print(out3)
cdf2=rbind(rep(0,3),out3)
cdf2=cbind(rep(0,4),cdf2)
pmf=apply(cdf2,2,diff)
pmf2=apply(t(pmf),2,diff)
pmf2=t(pmf2)  # rectangle probabilities
print(pmf2)

cat("\n pbvt rho changing\n")
z1=.3; z2=.4; rho=seq(-.9,.9,.1); nu=2
param=cbind(rho,rep(nu,length(rho)))
out1=pbvt(z1,z2,param)
print(cbind(rho,out1))
cat("\n pbvt z1 changing\n")
z1=seq(-2,2,.4)
z2=.4; rho=.5; nu=2
out2=pbvt(z1,z2,c(rho,nu))
print(cbind(z1,out2))
}

\keyword{distribution}
