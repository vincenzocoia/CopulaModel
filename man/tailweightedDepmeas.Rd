\name{tailweightedDepmeas}
\Rdversion{1.1}
\alias{tailweightedDepmeas}
\alias{twdm}
\alias{twdm.emp}
\alias{twdm.emp.vec}
\alias{twdmnestcop}
\alias{twdm1factcop}
\alias{twdm2factcop}
\title{
Tail-weighted dependence measures
}
\description{
Tail-weighted dependence measures: (a) bivariate copulas and (b) empirical
data
}
\usage{
twdm(pcop,param,power,nq,tscore=F)
twdm.emp(data,power)
twdm.emp.vec(data,power)  # vectorized, faster version
twdmnestcop(dcop,pcondcop,param1,param2,power,nq) # bivariate margin of
#                       nested-factor copula, 2 different groups
twdm1factcop(pcondcop,param,power,nq) # bivariate margin of 1-factor copula
twdm2factcop(pcondcop1,pcondcop2,param1,param2,power,nq) # bivariate margin of 2-factor copula
}
\arguments{
\item{pcop}{function for bivariate copula cdf}
\item{param}{dependence parameter of pcop, or pcondcop in the two variables}
\item{power}{power to use for tail-weighted dependence measure, good choice 
is 6 for twdm }
\item{nq}{number of quadrature points for Gauss-Legendre quadrature }
\item{tscore}{if T,
Student t transform or normal transform of Gauss-Legendre quadrature
points are used; this is faster is pcop is pbvncop or pbvtcop }
\item{data}{data matrix with dimensions nxd }
\item{dcop}{function for bivariate copula density for global with group1
and group2 latent}
\item{pcondcop}{function for copula conditional cdf given latent}
\item{pcondcop1}{function for copula conditional cdf for first latent}
\item{pcondcop2}{function for copula conditional cdf for second latent}
\item{param1}{dependence parameter of pcondcop1 in the two variables for
twdm2factcop; dependence parameter of dcop for two group variables for twdmnestcop}
\item{param2}{dependence parameter of pcondcop2 in the two variables;
dependence parameter of pcondcop for two observed variables in twdmnestcop}
}
\value{
\item{twdm}{for twdm, twdmnestcop, twdm1factcop, twdm2factcop: 
vector of length 2 lower and upper tail-weighted dependence measure values}
\item{ltwdm}{for twdm.emp and twdm.emp.vec: 
dxd matrix of empirical lower tail-weighted 
dependence measure values; for upper tail-weighted values, input with
negation of the data set.}
}
\seealso{
\code{\link{factorcopcdf}}
\code{\link{factorcopsim}}
\code{\link{pcond}}
\code{\link{structcop}}
}
\references{
Krupskii P (2014). Structured Factor Copulas and Tail Inference.
PhD thesis, University of British Columbia.

Krupskii P and Joe H (2015).
Tail-weighted measures of dependence.
J Applied Statistics, 42, 614-629.
}
\examples{
th1=gum.b2cpar(.7)
th2=gum.b2cpar(.6)
gum0.tw=twdm(pgum,th1,power=6,nq=15)
# 1-factor and 2-factor
gum1.tw=twdm(pfact1gum,c(th1,th1),power=6,nq=15)
gum2.tw=twdm(pfact2gum,matrix(c(th1,th1,th2,th2),2,2),power=6,nq=15)
gum1b=twdm1factcop(pcondgum,c(th1,th1),6,35)
gum2b=twdm2factcop(pcondgum,pcondgum,c(th1,th1),c(th2,th2),6,35)
# theoretical
cat(gum0.tw,"\n")
cat(gum1.tw,gum1b,"\n") # same from the two methods
cat(gum2.tw,gum2b,"\n") # same from the two methods
#
n=1000
set.seed(123)
gumdat1=sim1fact(n,c(th1,th1),qcondgum,"gumbel",ivect=TRUE)
set.seed(124)
gumdat2=sim2fact(n,c(th1,th1),c(th2,th2),qcondgum,qcondgum,"gumbel","gumbel",ivect=TRUE)
# empirical
gum1.ltw=twdm.emp(gumdat1,power=6)
gum1.utw=twdm.emp(1-gumdat1,power=6)
print(c(gum1.ltw[1,2],gum1.utw[1,2]))
gum2.ltw=twdm.emp.vec(gumdat2,power=6)
gum2.utw=twdm.emp.vec(1-gumdat2,power=6)
print(c(gum2.ltw[1,2],gum2.utw[1,2]))
# nested-factor
gumn=twdmnestcop(dgum,pcondgum,c(th1,th1),c(th2,th2),6,55) 
cat(gumn,"\n")
n=1000
set.seed(123)
grsize=c(2,2)
gumdatn=simnestfact(n,grsize, cop=3, c(th1,th1,th2,th2,th2,th2))
gumn.ltw=twdm.emp.vec(gumdatn,power=6)
gumn.utw=twdm.emp.vec(1-gumdatn,power=6)
print(c(gumn.ltw[1,2],gumn.utw[1,2]))
}
\keyword{copula}
\keyword{dependence}
