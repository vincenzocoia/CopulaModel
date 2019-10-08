\name{ordprobit.univar}
\alias{ordprobit.univar}
\alias{mord2uu}
\title{
Maximum likelihood for ordinal probit model
}
\description{
Maximum likelihood for ordinal probit: 
Newton-Raphson minimization of negative log-likelihood,
and conversion to uniform/normal scales for fitting copula model in case of
repeated measures with fixed cluster size
}
\usage{
ordprobit.univar(x,y,iprint=F,mxiter=20,toler=1.e-6)
mord2uu(xmat,yvec,nrep,b0cut,bvec) # multivariate ordinal to (0,1) vector
}
\arguments{
\item{x}{vector or matrix of explanatory variables. Each row corresponds
to an observation and each column to a variable.  The number of rows of x
should equal the number of data values in y, and there should be fewer 
columns than rows.  Missing values are not allowed.}
\item{y}{numeric vector containing the ordinal response. 
The values must be in the range 1,2,...,ncateg or 0,1,...(ncateg-1),
where ncateg is the number of categories.
Missing values are not allowed.}
\item{iprint}{print flag for the iterations 
for numerical maximum likelihood, default is FALSE}
\item{mxiter}{maximum number of Newton-Raphson iterations}
\item{toler}{tolerance for convergence in Newton-Raphson iterations}
\item{xmat}{vector or matrix of explanatory variables; like above x}
\item{yvec}{similar to above y}
\item{nrep}{number of repeated measures or cluster size for each subject/unit}
\item{b0cut}{vector of cutpoints}
\item{bvec}{vector of regression coefficients}
}
\value{For ordprobit.univar(),
list of MLE of parameters and their associated standard 
errors, in the order cutpt1,...,cutpt(number of categ-1),b1,...b(number of covariates).
$negloglik for value of negative log-likelihood, evaluated at MLE;
$cutpts for MLE of ordered cutpoint parameters;
$beta for MLE of regression parameters;
$cov for estimated covariance matrix of the parameters.

For mord2uu, a list with components
$uudat for transform of cdf to U(0,1),
$zzdat for transform of cdf to N(0,1).
}
\seealso{
\code{\link{ordinal}}
}
\details{
If ordprobit for repeated measures ordinal probit fails to converge
from the simple starting point in that function, this function
ordprobit.univar should provide a better starting point.
It is also equivalent to ordprobit with an identity latent
correlation matrix.

The ordinal probit model is similar to the ordinal logit model
(proportion odds logistic regression : polr in library MASS),
The parameter estimate of ordinal logit are roughly 1.8 to 2 times
those of ordinal probit (the signs of the parameters in polr
may be different, as this function may be using a different
orientation for the latent variable). 
}
\references{
Anderson JA and Pemberton JD (1985). The grouped continuous model for
multivariate ordered categorical variables and covariate adjustment.  
Biometrics, 41, 875-885.

}
\examples{
data(ordinalex)
xvec=c(t(ordinalex$xx))
yvec=c(t(ordinalex$yy))
ord.univar=ordprobit.univar(xvec,yvec,iprint=TRUE)
print(ord.univar)
ord.univar2=ordprobit.univar(xvec,yvec-1,iprint=TRUE)
print(ord.univar2) # same as ord.univar
ordtr=mord2uu(xvec,yvec,4,ord.univar$cutpts,ord.univar$beta)
ordtr2=mord2uu(xvec,yvec-1,4,ord.univar$cutpts,ord.univar$beta) #same
max(abs(ordtr$uudat-ordtr2$uudat))
}
\keyword{distribution}
\keyword{regression}
