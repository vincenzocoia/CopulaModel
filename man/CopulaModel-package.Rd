\name{CopulaModel-package}
\alias{CopulaModel-package}
\alias{CopulaModel}
\docType{package}
\title{
CopulaModel: Dependence Modeling with Copulas
}
\description{
CopulaModel is the accompanying software for the book:
Dependence Modeling with Copulas, by H. Joe, Chapman & Hall/CRC, 2014.
With this software, a reader can check (almost) all numerical
computations in the book. Much of the contributions for factor copula
models is also described in the Ph.D. thesis of P. Krupskii (2014).  
}
\details{
All of the algorithms in Dependence Modeling with Copulas are shown within
this software. There are templates for doing many numerical calculations
and maximum likelihood estimation with copulas. 
Included are examples of faster code with links to Fortran90
and C.  It is not possible to provide functions for all
possible uses of copulas. But a user can adapt the templates and code
within this software for many other applications.
The source code (R, Fortran90, C) has more documentation than these
help pages. The functions do not check on the domains of the arguments
such as the copula parameter, but the domains are indicated in
the source code. If a function doesn't seem to be working with your inputs,
check the source code and look for examples in the 
tests subdirectories and examples in the book chapter subdirectories. 

Some features include:
\itemize{
\item For the 1-, 2-, 3-parameter bivariate copula
families in the book, most have pcop, dcop, rcop, pcondcop,
qcondcop for the copula cdf, copula density, copula random
generation, conditional distribution \eqn{C_{2|1}} and 
conditional quantile \eqn{C_{2|1}^{-1}}, where 'cop' is an abbreviated
copula name.
Also for many bivariate copula families, there are
conversions among copula parameter,
Kendall's tau, Spearman's rho, Blomqvist's beta, correlation of normal scores,
and tail dependence parameters.

\item R-vine (regular vine) models for continuous data with specified vine 
array and parametric pair-copula families.
\item R-vine models for discrete response data with possibility of covariates.
\item Factor copula models for continuous response (factor models here
are truncated vines with latent variables).
\item Factor copula models for ordinal (item) response.
}
}
\author{
H. Joe and P. Krupskii. 
}
\references{
H. Joe (2014). Dependence Modeling with Copulas. Chapman & Hall/CRC.
Boca Raton, FL.

P. Krupskii (2014). Structured Factor Copulas and Tail Inference. 
PhD thesis, University of British Columbia.

Other specific references are given within code files or help pages.
}
\keyword{ package }
