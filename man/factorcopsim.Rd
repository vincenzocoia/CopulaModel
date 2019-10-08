\name{factorcopsim}
\Rdversion{1.1}
\alias{factorcopsim}
\alias{sim1fact}
\alias{sim2fact}
\title{
Simulation for factor copulas
}
\description{Simulation for 1-factor copula and 2-factor copula,
all bivariate linking copulas in same parametric family for each factor
}
\usage{
sim1fact(n,parobj1,qcond1,copname1,ivect=F)
sim2fact(n,parobj1,parobj2,qcond1,qcond2,copname1="a1",copname2="a2",ivect=F)
}
\arguments{
\item{n}{sample size }
\item{parobj1}{parameter vector of dimension d or parameter matrix with
d rows, where d is dimension of factor copula; parobj is dx2 for something
like BB1 copula }
\item{parobj2}{parameter vector of dimension d or parameter matrix with
d rows, for linking copulas to factor 2}
\item{qcond1}{name of conditional copula function \eqn{C_{U|V}^{-1}(u|v)},
choices include qfrk, qgum, qgumr, qbb1, qt with fixed nu1.  }
\item{qcond2}{name of conditional copula function \eqn{C_{U|V}^{-1}(u|v)}
for second factor,
choices include qfrk, qgum, qgumr, qt with fixed nu2. }
\item{copname1}{copula name: the function checks on "frank", "mtcj", "mtcjr", 
"fgm" for which qcond has closed form. }
\item{copname2}{copula name for factor 2}
\item{ivect}{flag that is T if qcond1 and qcond2 have vectorized forms}
}
\value{data matrix of dimension nxd
}
\seealso{
\code{\link{factorcopmle}}
\code{\link{rfactcop}}
}
\references{
Krupskii P and Joe H (2013).
Factor copula models for multivariate data.
Journal of Multivariate Analysis, 120, 85-101.
}
\examples{ 
bevec=c(.8,.7,.6,.5,.5)
cpar.frk=frk.b2cpar(bevec)
lmbb1=matrix(c(.3,.4,.5,.3,.5, .6,.6,.6,.7,.7),5,2)
cpar.bb1=lmbb1
for(i in 1:nrow(lmbb1))
{ cpar.bb1[i,]=bb1.lm2cpar(lmbb1[i,]) }
n=300
set.seed(123)
frkdat=sim1fact(n,cpar.frk,qcondfrk,"frank")
print(cor(frkdat))
set.seed(123)
bb1dat=sim1fact(n,cpar.bb1,qcondbb1,"bb1",ivect=FALSE)
print(cor(bb1dat))
set.seed(123)
bb1dat=sim1fact(n,cpar.bb1,qcondbb1,"bb1",ivect=TRUE)
print(cor(bb1dat))
# pairs(bb1dat)
bevec2=c(.0001,.6,.6,.6,.7)
cpar.frk2=frk.b2cpar(bevec2)
frk2dat=sim2fact(n,cpar.frk,cpar.frk2,qcond1=qcondfrk,qcond2=qcondfrk,"frank","frank")
print(cor(frk2dat))
}
\keyword{copula}
\keyword{factor model}
\keyword{simulation}
