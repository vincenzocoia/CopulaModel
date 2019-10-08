\name{extremevalue}
\Rdversion{1.1}
\alias{gevmle}
\alias{gpmle}
\alias{dgev}
\alias{pgev}
\alias{qgev}
\alias{logdgev}
\alias{dgpareto}
\alias{pgpareto}
\alias{qgpareto}
\title{
Maximum likelihood for GEV and generalized Pareto
}
\description{
Maximum likelihood for generaluzed extreme value and generalized Pareto, using
modified Newton-Raphson
}
\usage{
gevmle(xdata, maxitn=20)
gpmle(xdata ,maxitn=20)
dgev(x, xi=1,mu=0,sigma=1) 
pgev(x, xi=1,mu=0,sigma=1)
qgev(p, xi=1,mu=0,sigma=1)
logdgev(x, xi=1,mu=0,sigma=1)
dgpareto(x,xi, sigma=1)
pgpareto(x,xi, sigma=1)
qgpareto(p,xi, sigma=1)
}
\arguments{
\item{xdata}{data set, should be positive-valued exceedances for gpmle}
\item{maxitn}{maximum number of iterations for Newton-Raphson method}
\item{x}{scalar of vector}
\item{p}{scalar of vector, values in 0 to 1}
\item{xi}{tail index or shape parameter of GEV}
\item{mu}{location parameter}
\item{sigma}{scale parameter}
}
\value{list with $loglik, $params,  $covar for gevmle and gpmle functions

density, cdf, quantile or log density values with the other functions
}
\examples{
set.seed(123)
x=rnorm(2000)
xmat=matrix(x,40,50)
mxdat=apply(xmat,1,max)
gevmle(mxdat) # xi is negative because of sub-asymptotics
gpmle(x[x>1.3]-1.3)
rpareto=function(n,alp,s)
{ u=runif(n); tem=u^(-1/alp)-1; return(s*tem) }
set.seed(123)
x=rpareto(2000,alp=3,s=1)
xmat=matrix(x,40,50)
mxdat=apply(xmat,1,max)
gevmle(mxdat) # xi close to 1/alp=1/3
gpmle(x[x>1.2]-1.2)
}
\keyword{distribution}
\keyword{maximum likelihood}
