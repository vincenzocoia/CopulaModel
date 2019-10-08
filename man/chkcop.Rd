\name{chkcop}
\Rdversion{1.1}
\alias{chkcopderiv}
\alias{chkcopcond}
\alias{chkincrease}
\title{
Checks for pcop pcond dcop qcond functions
}
\description{
Check pcop pcond dcop qcond functions for copula cdf, conditional cdf,
pdf, conditional quantile
}
\usage{
chkcopderiv(u,vvec,cpar,bcdf,pcond,bpdf,str=" ",eps=1.e-4)
chkcopcond(uvec,vvec,cpar,pcond,qcond,str=" ",tol=1.e-5)
chkincrease(mat)
}
\arguments{
\item{u}{scalar in (0,1)}
\item{vvec}{vector of values in (0,1)}
\item{uvec}{vector of values in (0,1), same length as vvec}
\item{cpar}{copula parameter}
\item{bcdf}{function for bivariate copula cdf = pcop}
\item{bpdf}{function for bivariate copula pdf = dcop}
\item{pcond}{function for conditional copula cdf}
\item{qcond}{function for conditional copula quantile}
\item{str}{string for copula name}
\item{eps}{step size for numerical derivative in chkcopderiv}
\item{tol}{tolerance for pcond/qcond composition versus identity in chkcopcond}
\item{mat}{matrix or dataframe}
}
\value{
For chkcopderiv and chkcopcond,
output is printed, nothing is returned.
chkcopcond checks pcond vs numerical derivative of pcop;
chkcopderiv checks also dcop vs numerical derivative of pcond.

For chkincrease, a value of 1 means the values are increasing in
each row, otherwise 0 is returned.
}
\examples{
# checks cdf, conditional cdf, pdf for Plackett copula
be=.8
par.pla=pla.b2cpar(be)
u=.3
vvec=seq(.4,.9,.1)
chkcopderiv(u,vvec,par.pla,bcdf=ppla,pcond=pcondpla,bpdf=dpla,str="pla")
uvec=seq(.1,.6,.1)
vvec=seq(.4,.9,.1)
chkcopcond(uvec,vvec,par.pla,pcondpla,qcondpla,"pla")
# check monotonicity for Plackett copula
u=seq(.1,.9,.2)
umat=matrix(rep(u,5),5,5,byrow=TRUE)
vmat=matrix(rep(u,5),5,5)
cdf=ppla(umat,vmat,par.pla)
chkincrease(cdf)
# check monotonicity of dependence measure for Galambos and Plackett copulas
data(deptabder)
chkincrease(gal.deptab[,1:6]) # cols cpar beta tau rhoS rhoN lambda
apply(gal.deptab,2,diff)
chkincrease(pla.deptab[,1:5]) # cols cpar beta tau rhoS rhoN l
}
\keyword{copula}
