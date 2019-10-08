\name{rcop}
\Rdversion{1.1}
\alias{rcop}
\alias{rpla}
\alias{rfrk}
\alias{rmtcj}
\alias{rjoe}
\alias{rgum}
\alias{rgumr}
\alias{rgal}
\alias{rhr}
\alias{rbreflasym}
\alias{rbpermasym}
\alias{rbasymgum1}
\alias{rbMO1}
\title{
Simulation from parametric bivariate copula families
}
\description{
Simulation from parametric bivariate copula families (one-parameter and
two-parameter)
}
\usage{
#rcop(n,cpar)  
rcop(n,cpar,icheck=F) 
}
\arguments{
\item{n}{sample size}
\item{cpar}{copula parameter: could be scalar or vector depending on the 
copula family}
\item{icheck}{if T, output average u, average v and (Spearman) correlation}
}
\details{
Choices are 'cop' in rcop 

pla, frk, mtcj, joe, gum, gumr (reflected Gumbel), gal, hr, 
bb1, bb2, bb3, bb4, bb5, bb6, bb7, bb8, bb9, bb10 for bivariate families;

See help page for pcop for the abbreviations of the copula names.

breflasym, bpermasym, basymgum1, bMO1 for 'cop' in rcop
for extreme cases of bivariate reflection or permutation asymmetric:

bivariate permutation asymmetric,

bivariate asymmetric Gummbel with one skew parameter,

bivariate Marshall-Olkin skew 1-parameter subfamily.
}
\value{
nx2 matrix with values in (0,1)
}
\seealso{
\code{\link{pcop}}
\code{\link{rbivcop2param}}
}
\examples{
n=1000
set.seed(12345)
cpar=gum.tau2cpar(.45)
udat=rgum(n,cpar)
taucor(udat[,1],udat[,2])
cpar=depmeas2cpar(.45,"rhoS","gumbel")
udat=rgum(n,cpar)
cor(udat[,1],udat[,2],method="spearman")
cor(udat[,1],udat[,2])
# see help page for skewrefl for use of 
# rbreflasym, rbpermasym, rbsymgum1, rbMO1
}
\keyword{copula}
\keyword{simulation}
