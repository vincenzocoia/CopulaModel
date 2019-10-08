\name{negbinom}
\Rdversion{1.1}
\alias{negbinom}
\alias{nb1nllk}
\alias{nb2nllk}
\alias{nbpmfcdf}
\alias{nb1pmfcdf}
\alias{nb2pmfcdf}
\alias{nb1cdf}
\alias{nb2cdf}
\alias{nb1pmf}
\alias{nb2pmf}
\title{
Negative binomial regression for count data
}
\description{
Negative binomial (NB) regression for count data; 2 versions.
NB(theta,xi) has pmf 
f(y;theta,xi)= [ Gamma(theta+y) xi^y] / [ Gamma(theta) y! (1+xi)^(theta+y) ]
}
\usage{
nb1nllk(param,y,xdat)
nb2nllk(param,y,xdat)
nbpmfcdf(ub,theta,p)
nb1pmfcdf(ub,param,x)
nb2pmfcdf(ub,param,x)
nb1cdf(y,param,x)
nb2cdf(y,param,x)
nb1pmf(y,param,x)
nb2pmf(y,param,x)
}
\arguments{
\item{param}{
parameter of NB model, length is 2+number of covariates; the parameters are: 
b0=intercept, bvec= vector regression coefficients
(length(bvec)=length(x)=ncol(xdat), and finally xi or theta.
For NB1, mu(x)= exp(b[0]+bvec^T x), xi=(overdispersion index minus
one) is fixed, and theta(x)=mu(x)/xi.
For NB2, mu(x)= exp(b[0]+bvec^T x), theta=convolution parameter
is fixed and xi(x)=mu(x)/theta.}
\item{theta}{convolution parameter}
\item{p}{probability parameter between 0 and 1}
\item{xdat}{matrix for nb1nllk and nb2nllk} 
\item{x}{vector for nb1pmfcdf, nb2pmfcdf, nb1cdf, nb2cdf, nb1pmf, nb2pmf}
\item{y}{vector for nb1nllk and nb2nllk (with length(y)=nrow(xdat));
non-negative integer for the other functions}
\item{ub}{upper bound integer for which pmf and cdf are computed}
}
\value{
negative log-likelihood  for nb1nllk and nb2nllk;
matrix with columns (0:ub,pmf,cdf) for nbpmfcdf, nb1pmfcdf and nb2pmfcdf,
computed in an efficient way (parameters assumed to be such that
most probability is on small counts);
cdf for nb1cdf and nb2cdf;
pmf for nb1pmf and nb2pmf.
}
\seealso{
\code{\link{gpoisson}}
}

\examples{
y= c(
 2, 1, 1, 0,35, 9, 0, 1, 4, 0, 0, 1, 4, 0, 0, 8, 7, 2, 0, 7, 0, 0, 3, 4, 0,
 4, 1, 3, 0, 6, 1, 0, 2, 8, 0,12, 0, 4, 2, 1, 3, 0, 9, 0, 0, 0, 2, 0, 8, 1,
 2, 4, 2, 0, 0, 2, 1, 3, 2, 1, 3, 4, 4, 5, 0, 4, 0, 2, 0,28, 1,24, 1, 0,10,
 3, 3, 0, 0, 7, 2, 4, 6, 4,13, 5, 8, 0, 1, 6, 0,24, 9, 0,10, 0, 0, 8, 5, 3,
16, 0, 4, 1, 1, 4,12, 4, 3, 5, 0, 2, 1, 5, 3, 0, 0, 6, 4, 2, 0, 2, 0,15, 3,
 0, 2, 3, 4, 5, 0, 3, 0, 0, 6, 0, 0,15, 0, 0, 0, 1, 3, 0, 1, 0, 4, 2,10, 4,
 1, 0, 0, 0, 5, 0, 0, 2, 0, 4, 0, 0, 2,25, 0, 0,13, 0, 0,21, 3, 0, 0, 0, 2,
 2, 0, 4,13, 2, 9, 9, 2, 0, 1, 2, 2, 8, 6, 0, 4, 1, 2, 0, 0, 0, 0, 0, 0, 2,
 2, 0, 3, 1, 1, 7, 3, 0, 2, 2, 1, 3, 2, 2, 1, 3, 3, 0, 0, 0, 2, 0, 0, 0, 0,
 1, 2, 2, 0, 0, 9, 0, 0, 1, 1, 0, 2,10, 0,17, 2, 0,14, 0, 5, 9, 2, 0, 6, 3,
 3, 1, 0,11, 4, 9, 0, 1, 0, 0,12, 4, 0, 1,21, 0, 3, 2, 0, 1, 0, 1, 3, 8,10,
19, 0, 2, 7, 1, 0, 2, 0, 4, 0, 6, 4, 7, 1, 0, 1, 3, 4, 0, 4)
hsat=c(
 8, 7, 3,10, 6, 5, 8, 9, 9, 8,10, 8, 6, 7,10, 8, 5, 8, 8, 6, 8, 8, 8, 9,10,
 7, 9,10, 8, 6, 6, 9, 7, 5,10, 4, 8, 4, 5, 5, 7, 6, 7,10, 9, 9, 5, 7, 4, 7,
 6, 6, 7, 5,10, 9,10, 7, 8, 6, 5, 5, 0, 5, 7, 3, 8, 8, 7, 5, 5, 0, 7, 6, 3,
10, 7, 7,10, 5, 5, 4, 2, 7, 6, 2, 5,10, 7, 8, 5, 5, 5,10, 3, 9, 6, 8,10,10,
 4, 7, 2, 8, 9, 0, 0, 5, 8, 3, 7, 6,10, 4, 5, 7, 6, 7, 3, 4,10, 4, 8, 8, 3,
 9, 5,10, 9, 5,10,10, 8,10, 5,10, 6, 5, 9, 8,10, 7, 8, 9, 7, 8, 4, 8, 3, 5,
 5, 7,10, 8, 1, 3, 3, 8,10, 3, 5, 5, 7, 5,10, 8, 5, 8, 5, 0, 6, 8, 2, 5, 6,
 7,10, 5, 0, 5, 2, 0, 3,10, 7, 4, 6, 9, 2, 8, 5, 9, 7, 5,10, 8, 8, 7, 7, 7,
10,10, 2, 5, 7, 5, 9, 6, 7, 6, 9, 9, 6, 8,10, 7, 8, 8,10,10, 5,10, 5, 8,10,
 8, 7,10, 9,10, 4, 6, 9, 5, 9, 9, 6, 8, 8, 2, 5, 8, 3, 7, 0, 8, 8,10, 5, 7,
 6, 7,10, 5, 5, 1, 5, 6, 4,10, 5, 5, 5, 7, 2, 8, 5,10,10,10,10, 6, 6, 6, 6,
 7, 8, 8,10,10, 8, 7, 8, 3, 8, 8, 8, 6, 3, 7,10,10, 2, 9, 2)
fit1=nlm(ieenllk,p=c(2.5,-.2,4),hessian=TRUE,print.level=1,upmf=nb1pmf,
    xdat=hsat,ydat=y,LB=c(-1,-2,0),UB=c(10,10,10))
fit2=nlm(ieenllk,p=c(2.5,-.2,0.8),hessian=TRUE,print.level=1,upmf=nb2pmf,
    xdat=hsat,ydat=y,LB=c(-1,-2,0),UB=c(10,10,10))
}
\keyword{distribution}
\keyword{regression}
