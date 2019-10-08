\name{asymgum}
\Rdversion{1.1}
\alias{asymgum}
\alias{Basymgum}
\alias{pasymgum}
\alias{dasymgum}
\alias{pcondasymgum12}
\alias{pcondasymgum21}
\alias{asymgum.cpar2tau}
\alias{asymgum.cpar2rhoS}
\alias{asymgal}
\alias{Basymgal}
\alias{pasymgal}
\alias{dasymgal}
\alias{pcondasymgal12}
\alias{pcondasymgal21}
\alias{asymgal.cpar2tau}
\alias{asymgal.cpar2rhoS}
\alias{asymgumMO}
\alias{AasymgumMO}
\alias{pasymgumMO}
\alias{dasymgumMO}
\alias{pcondasymgumMO12}
\alias{pcondasymgumMO21}
\alias{asymgumMO.cpar2tau}
\alias{asymgumMO.cpar2rhoS}
\title{
Bivariate asymmetric Gumbel and Galambos copulas
}
\description{
Bivariate asymmetric Gumbel and Galambos copulas
}
\usage{
pasymgum(u,v,cpar)  # C(u,v)
pcondasymgum21(v,u,cpar) # C_{2|1}(v|u)
pcondasymgum12(u,v,cpar)  # C_{1|2}(u|v))
dasymgum(u,v,cpar)  # c(u,v)
asymgum.cpar2tau(cpar) # Kendall's tau
asymgum.cpar2rhoS(cpar) # Spearman's rho
Basymgum(ww,cpar,mxiter=30,eps=1.e-7,iprint=F)  # B(w)=A(w,1-w)
# similar to above with 'gum' replaced by 'gal' or 'gumMO' except
AasymgumMO(x,y,cpar)
}
\arguments{
\item{u}{value in interval 0,1; could be a vector }
\item{v}{value in interval 0,1; could be a vector }
\item{cpar}{copula parameter vector, 
for asymgum, the two parameters are each in (0,1);
for asymgal, the two parameters are each negative;
for asymgumMO, the 3 parameters are delta>1, pi1 in (0,1] and pi2 in (0,1]}
\item{ww}{vector of values in (0,1)}
\item{mxiter}{maximum number of iterations}
\item{eps}{tolerance for convergence}
\item{iprint}{print flag for intermediate results}
\item{x}{vector of positive values}
\item{y}{vector of positive values}
}
\value{
cdf for pasymgum

conditional cdf for pcondasymgum21 and pcondasymgum12

pdf for dasymgum

Kendall's tau for asymgum.cpar2tau

Spearman's rho for asymgum.cpar2rhoS

list with $Bfn (function values), $Bder (first derivatives), 
$Bder2 (second derivatives) for Basymgum
}
\details{
The distribution in asymgum is also known as bilogistic in the multivariate 
extreme value literature -- it is the bivariate Gumbel when the
two parameters are equal and the Gumbel parameter is the reciprocal.
The distribution in asymgal is also known as negative bilogistic in 
the multivariate extreme value literature
-- it is the bivariate Galambos when the
two parameters are equal and the Galambos parameter is the negative reciprocal.
The distribution in asymgumMO becomes bivariate Gumbel when the second and
third parameters are both 1.
}
\examples{
ze=.3; eta=.2
cpar=c(ze,eta)
heps=1.e-5
u=.3; v=.8
chkcopderiv(u,v,cpar,bcdf=pasymgum,pcond=pcondasymgum21,bpdf=dasymgum,
  str="asymgum",eps=heps)
tau=asymgum.cpar2tau(cpar)
rho=asymgum.cpar2rhoS(cpar)
cat(cpar,tau,rho,"\n")
# special case of Gumbel
cpar=c(.5,.5)
cdf=pasymgum(.4,seq(.1,.9,.2),cpar)
cdf2=pgum(.4,seq(.1,.9,.2),1/cpar[1])
print(cbind(cdf,cdf2))
tau=asymgum.cpar2tau(cpar)
rho=asymgum.cpar2rhoS(cpar)
cat(tau,rho,"\n")
tau=gum.cpar2tau(1/cpar[1])
rho=gum.cpar2rhoS(1/cpar[1])
cat(tau,rho,"\n")
# asymmetric Galambos
cpar=c(-.5,-.5)
cdf=pasymgal(.4,seq(.1,.9,.2),cpar)
cdf2=pgal(.4,seq(.1,.9,.2),-1/cpar[1])
print(cbind(cdf,cdf2))
#and gumMO to add
}

\keyword{copula}
