\name{KLdiv}
\Rdversion{1.1}
\alias{KLcopvscop}
\alias{KLcopvsbvn}
\alias{KL12gl}
\alias{KLoptgl}
\title{
Kullback-Leibler divergence and sample size for two bivariate copula
densities
}
\description{
Kullback-Leibler divergence and sample size for two bivariate copula
densities
}
\usage{
KLcopvsbvn(rho,dcop2,param2,copname2="bivcop",UB=7,iprint=F)
KLcopvscop(copname1="cop1",param1,dcop1,copname2="cop2",param2,dcop2,UB=7,iprint=F)
KL12gl(par2,dcop2,par1,dcop1,parlb=0,gl)
KLoptgl(dcop1,dcop2,par1,par2,name1,name2,par2lb=0,gl,pcop1,pcop2,
   ccdf1a,ccdf1b,ccdf2a,ccdf2b,prlevel=0)
}
\arguments{
\item{dcop1}{function for first bivariate copula density; form is dcop(u,v,param)}
\item{dcop2}{function for second bivariate copula density}
\item{copname1}{name of first bivariate copula density}
\item{copname2}{name of second bivariate copula density}
\item{rho}{parameter of dcop1 if it is bivariate normal}
\item{param1}{parameter of dcop1}
\item{param2}{parameter of dcop2}
\item{UB}{limit to use for 2-dimensional integration with respect to 
bivariate normal density}
\item{iprint}{print flag for intermediate results, FALSE by default}
\item{par1}{copula parameter for family 1}
\item{par2}{copula parameter for family 2}
\item{parlb}{lower bound for copula parameter for family 2}
\item{par2lb}{lower bound for copula parameter for family 2}
\item{gl}{Gauss-Legendre quadrature object with $nodes and $weights}
\item{name1}{name of first bivariate copula density}
\item{name2}{name of second bivariate copula density}
\item{pcop1}{function for first bivariate copula cdf; form is pcop(u,v,param)}
\item{pcop2}{function for second bivariate copula cdf}
\item{ccdf1a}{function for \eqn{C_{2|1}} for pcop1}
\item{ccdf1b}{function for \eqn{C_{1|2}} for pcop1}
\item{ccdf2a}{function for \eqn{C_{2|1}} for pcop2}
\item{ccdf2b}{function for \eqn{C_{1|2}} for pcop2}
\item{prlevel}{print.level for nlm}
}
\value{For KLcopvsbvn and KLcopvscop, a vector with 5 elements:
KLcop1true, KLcop2true, Jeffreys, sampsize1, sampsize2,
where KLcop1true (KLcop2true) is the KL divergence when
copula 1 (respectively copula 2) is the true model;
Jeffreys=KLcop1true+KLcop2true is Jeffreys' divergence,
sampsize1 (sampsize2) is the sample size needed to distinguish with
the other copula family wth probability 0.95 when 
copula 1 (respectively copula 2) is the true model.

For KL12gl, a single value for divergence; this function is optimized
by KLoptgl.

For KLoptgl, a list with
$cpar2 for the parameter of copula2 leading to small KL divergence
with copula1 with parameter param1,
$depm1=c(be1,tau1,rhoS1,rhoN1) is a vector of beta,tau,rhoS,rhoN for
copula1 with parameter param1,
$depm2=c(be2,tau2,rhoS2,rhoN12 is a vector of beta,tau,rhoS,rhoN for
copula2 with parameter cpar2.
}
\seealso{
\code{\link{pcop}}
\code{\link{pcond}}
\code{\link{deppar2taurhobetalambda}}
}
\examples{
rho=bvn.b2cpar(0.5)
par.gum=gum.b2cpar(0.5)
par.gal=gal.b2cpar(0.5)
KLcopvsbvn(rho,dcop2=dgum,param2=par.gum,copname2="Gumbel",UB=6,iprint=TRUE) 
KLcopvscop("Gumbel",par.gum,dgum,"Galambos",par.gal,dgal,UB=6,iprint=TRUE)
# parameter from Plackett family closed to bivariate Gaussian copula with 
#   rhoS=0.5 
cpar1=bvn.rhoS2cpar(.5)
cpar2=pla.rhoS2cpar(.5)  # 5.115661
gl25=gausslegendre(25)
kl=KLoptgl(dbvncop,dpla,cpar1,cpar2,"bvn","plackett",0,gl25,pbvncop,ppla,
   pcondbvncop,pcondbvncop,pcondpla,pcondpla)
print(kl$cpar2)  # 4.671857
}
\keyword{copula}
