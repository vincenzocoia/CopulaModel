\name{exchmvn}
\alias{exchmvn}
\alias{exchmvn.deriv.margin}
\alias{exchmvn.deriv.rho}
\alias{pmnorm}
\title{
Exchangeable (positive) multivariate normal
}
\description{
Rectangle probability and derivatives of positive exchangeable multivariate
normal, and trivariate normal}
\usage{
exchmvn(lb,ub,rho, mu=0,scale=1,eps=1.e-06)
exchmvn.deriv.margin(lb,ub,rho,k,ksign, eps=1.e-06)
exchmvn.deriv.rho(lb,ub,rho, eps=1.e-06)
pmnorm(lb,ub,mu,sigma, eps=1.e-05)
}
\arguments{
\item{lb}{vector of lower limits of integral/probability  }
\item{ub}{vector of upper limits of integral/probability }
\item{rho}{correlation (positive constant over pairs)}
\item{mu}{mean vector }
\item{scale}{standard deviation }
\item{eps}{tolerance for numerical integration }
\item{k}{ margin for which derivative is to be taken, that is,
derivative of exchmvn(lb,ub,rho) with respect to lb[k] or ub[k];
use exchmvn.deriv.rho for derivative of exchmvn(lb,ub,rho) with respect to rho}
\item{ksign}{
  value is -1 for derivative of exchmvn(lb,ub,rho) with respect to lb[k],
  value is +1 for derivative of exchmvn(lb,ub,rho) with respect to ub[k] }
\item{sigma}{covariance matrix }
}
\value{
rectangle probability or a derivative
}
\details{The positive exchangeable multivariate normal distribution
has a stochastic representation as a one-factor model from which
rectangle probabilities can be written as 1-dimensional integrals.
pmnorm() from  Schervish (1984) is recommended only for dimension d=3; 
otherwise use pmvnorm() in library mvtnorm.
}
\references{
Kotz S and Johnson NL (1972). Continuous Multivariate Distributions. Wiley, 
New York, page 48.

Schervish M (1984). Multivariate normal probabilities with error bound.
Applied Statistics, 33, 81-94.
}
\seealso{
\code{\link{fact1mvn}}
}

\examples{
# The tests here show clearly what the function parameters are.
# step size for numerical derivatives (accuracy of exchmvn etc about 1.e-6)
heps = 1.e-4

cat("case 1: m=3\n")
m=3
a=c(-1,-1,-1)
b=c(2,1.5,1)
rho=.6
pr=exchmvn(a,b,rho)
cat("pr=exchmvn(avec,bvec,rho)=",pr,"\n")
cat("derivative wrt rho\n")
rho2=rho+heps
pr2=exchmvn(a,b,rho2)
drh.numerical= (pr2-pr)/heps
drh.analytic= exchmvn.deriv.rho(a,b,rho)
cat("   numerical: ", drh.numerical, ", analytic: ", drh.analytic,"\n")

cat("derivative wrt a_k,b_k, k=1,...,",m,"\n")
for(k in 1:m)
{ cat("  k=", k, " lower\n")
  a2=a
  a2[k]=a[k]+heps
  pr2=exchmvn(a2,b,rho)
  da.numerical = (pr2-pr)/heps 
  da.analytic= exchmvn.deriv.margin(a,b,rho,k,-1)
  cat("   numerical: ", da.numerical, ", analytic: ", da.analytic,"\n")
  cat("  k=", k, " upper\n")
  b2=b
  b2[k]=b[k]+heps
  pr2=exchmvn(a,b2,rho)
  db.numerical = (pr2-pr)/heps 
  db.analytic= exchmvn.deriv.margin(a,b,rho,k,1)
  cat("   numerical: ", db.numerical, ", analytic: ", db.analytic,"\n")
}


cat("\ncase 2: m=5\n")
m=5
a=rep(-1,m)
b=c(2,1.5,1,1.5,2)
rho=.6
pr=exchmvn(a,b,rho)
cat("pr=exchmvn(avec,bvec,rho)=",pr,"\n")
cat("derivative wrt rho\n")
rho2=rho+heps
pr2=exchmvn(a,b,rho2)
drh.numerical= (pr2-pr)/heps
drh.analytic= exchmvn.deriv.rho(a,b,rho)
cat("   numerical: ", drh.numerical, ", analytic: ", drh.analytic,"\n")

cat("derivative wrt a_k,b_k, k=1,...,",m,"\n")
for(k in 1:m)
{ cat("  k=", k, " lower\n")
  a2=a
  a2[k]=a[k]+heps
  pr2=exchmvn(a2,b,rho)
  da.numerical = (pr2-pr)/heps 
  da.analytic= exchmvn.deriv.margin(a,b,rho,k,-1)
  cat("   numerical: ", da.numerical, ", analytic: ", da.analytic,"\n")
  cat("  k=", k, " upper\n")
  b2=b
  b2[k]=b[k]+heps
  pr2=exchmvn(a,b2,rho)
  db.numerical = (pr2-pr)/heps 
  db.analytic= exchmvn.deriv.margin(a,b,rho,k,1)
  cat("   numerical: ", db.numerical, ", analytic: ", db.analytic,"\n")
}

}
\keyword{distribution}
