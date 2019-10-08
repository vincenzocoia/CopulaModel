\name{deptab}
\alias{deptabder}
\alias{deptab}
\alias{pla.deptab}
\alias{frk.deptab}
\alias{mtcj.deptab}
\alias{joe.deptab}
\alias{gum.deptab}
\alias{gal.deptab}
\alias{hr.deptab}
\title{Tables of dependence measures for 1-parameter bivariate copula
families}
\usage{
data(deptabder)
}
\description{
Tables (dataframes) in deptabder are
pla.deptab, frk.deptab, mtcj.deptab, joe.deptab, gum.deptab,
gal.deptab, hr.deptab with columns:
cpar, beta, tau, rhoS, rhoN, lambda, tauder, rhoSder, rhoNder.

The bivariate copula families are Plackett, Frank,
MTCJ=Mardia-Takahasi-Cook-Johnson,
Joe-B5, Gumbel, Galambos, Huesler-Reiss.

cpar=copula parameter,
beta=Blomqvist's beta,
tau=Kendall's tau,
rhoS=Spearman's rho,
rhoN=correlation of normal scores,
lambda=tail dependence parameter,
tauder=(estimated) derivative of tau with respect to beta,
rhosder=(estimated) derivative of rhos with respect to beta,
rhoNder=(estimated) derivative of rhoN with respect to beta.
The latter three are used for pcinterpolate() for interpolation in the
function depmeas2cpar (dependence measure to copula parameter).
Given one of tau, rhoS or rhoN, the corresponding value of beta
is interpolated from these tables, and then inverted to get the
copula parameter. 
}
\format{
Tables included are the following.
\describe{
\item{\code{pla.deptab}}{Dependence measure table for bivariate Plackett copula}
\item{\code{frk.deptab}}{Dependence measure table for bivariate Frank copula}
\item{\code{mtcj.deptab}}{Dependence measure table for bivariate MTCJ copula}
\item{\code{joe.deptab}}{Dependence measure table for bivariate Joe copula}
\item{\code{gum.deptab}}{Dependence measure table for bivariate Gumbel copula}
\item{\code{gal.deptab}}{Dependence measure table for bivariate Galambos copula}
\item{\code{hr.deptab}}{Dependence measure table for bivariate Huesler-Reiss copula}
}
}

\seealso{
\code{\link{pcinterpolate}}
\code{\link{depmeas2cpar}}
}
\keyword{dependence}
