\name{genbeta2}
\Rdversion{1.1}
\alias{genbeta2}
\alias{dbgb2}
\alias{pbgb2}
\alias{dbgb2cop}
\alias{pbgb2cop}
\alias{pcondbgb2cop}
\alias{dmgb2}
\alias{dmgb2cop}
\alias{logdmgb2cop}
\alias{rmgb2cop}
\alias{bgb2.cpar2lm}
\alias{mgb2.cpar2cor}
\title{
Copula based on generalized beta of order 2
}
\description{
Copula based on generalized beta of order 2
}
\usage{
dbgb2(y1,y2,param) # bivariate density
pbgb2(y1,y2,param)  # bivariate cdf
dbgb2cop(u,v,param) # bivariate copula density
pbgb2cop(u,v,param)  # bivariate copula cdf
pcondbgb2cop(v,u,param)  # bivariate copula conditional cdf
dmgb2(yvec,param)    # multivariate density
dmgb2cop(uvec,param)    # multivariate copula density
logdmgb2cop(uvec,param) # log of multivariate copula density
rmgb2cop(n,param)  # simulation
bgb2.cpar2lm(param) # upper tail dependence 
mgb2.cpar2cor(param) # correlation matrix if zeta>2
}
\arguments{
\item{y1}{positive value; could be a vector }
\item{y2}{positive value; could be a vector }
\item{yvec}{positive vector of length d}
\item{u}{value in interval (0,1); could be a vector }
\item{v}{value in interval (0,1); could be a vector }
\item{uvec}{vector of length d; values in (0,1)}
\item{param}{parameter vector of dimension d+1; (eta1,eta2,zeta) for bivariate 
and (etavec,zeta) for d-variate; all parameters are positive }
\item{n}{simulation sample size}
}
\value{
pdf for dbgb2, dmgb2, dbgb2cop, dmgb2cop;

cdf for pbgb2, pbgb2cop;

conditional cdf for pcondbgb2cop;

log density for logdmgb2cop (use for maximum likelihood);

random d-vectors for rmgb2cop;

upper tail dependence parameter for bgb2.cpar2lm;

correlation matrix for mgb2.cpar2cor
}
\references{
Yang X, Frees EW, Zhang Z (2011). A generalized beta copula with
applications in modeling multivariate long-tailed data.
Insurance: Mathematics and Economics 49, 265-284.
}

\examples{
y1=2; y2=3;
param=c(1.5,1.8,2.3)
dbgb2(y1,y2,param)
dmgb2(c(y1,y2),param) # same as above
u1=.2; u2=.3;
dbgb2cop(u1,u2,param)
dmgb2cop(c(u1,u2),param) # same as above
n=1000
param=c(4.5,4.8,2.3)
set.seed(124)
udat=rmgb2cop(n,param)
print(cor(udat))
cat("compare theoretical\n")
print(mgb2.cpar2cor(param))
dmgb2copnllk= function(param,udat)
{ n=nrow(udat)
  if(any(param<=0)) return(1.e10)
  nllk=0
  for(i in 1:n) nllk=nllk-logdmgb2cop(udat[i,],param)
  nllk
}
out=nlm(dmgb2copnllk,p=param,udat=udat,hessian=TRUE,print.level=1,iterlim=20)
print(sqrt(diag(solve(out$hess))))
\dontrun{
# contour of density with N(0,1) margins
zvec=seq(-3,3,.2)
contourBivCop(param,zvec,dcop=dbgb2cop) }
}
\keyword{copula}
\keyword{distribution}
