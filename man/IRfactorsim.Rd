\name{IRfactorsim}
\Rdversion{1.1}
\alias{IRfactorsim}
\alias{sim1irfact}
\alias{sim2irfact}
\title{
Simulation for factor copulas for ordinal item response data
}
\description{
Simulation for factor copulas for ordinal item response data,
all bivariate linking copulas in same parametric family for each factor
}
\usage{
sim1irfact(ucuts,n,parobj1,qcond1,copname1,ivect=F)
sim2irfact(ucuts,n,parobj1,parobj2,qcond1,qcond2,copname1="",copname2="",ivect=F)
}
\arguments{
\item{ucuts}{(ncat-1)xd matrix of cutpoints, 
increasing in each column and bounded in (0,1); ncat=number of categories}
\item{n}{sample size}
\item{parobj1}{parameter vector of dimension d or parameter matrix with
d rows for factor 1, where d is dimension of factor copula; this is dx2 for 
something like BB1 copula}
\item{parobj2}{parameter vector of dimension d or parameter matrix with
d rows for factor 2}
\item{qcond1}{function for copula conditional inverse cdf 
\eqn{C_{U|V}^{-1}(u|v)},
choices include qcondfrk, qcondgum, qcondgumr, qcondbb1, qcondbvtcop with fixed nu1.}
\item{qcond2}{function for copula conditional inverse cdf 
\eqn{C_{U|V}^{-1}(u|v)}
for second factor,
choices include qcondfrk, qcondgum, qcondgumr, qcondbvtcop with fixed nu2. }
\item{copname1}{copula name: the function checks on "frank", "mtcj", "mtcjr", 
"fgm" for which qcond has closed form. }
\item{copname2}{copula name for factor 2}
\item{ivect}{flag that is T if qcond1 and qcond2 have vectorized forms}
}
\value{data matrix of dimension nxd with values in 0,...,ncat-1
}
\seealso{
\code{\link{IRfactormle}}
\code{\link{factorcopsim}}
}
\references{
Nikoloulopoulos A K and Joe H (2015).
Factor copula models for item response data, Psychometrika, 80, 126-150.
}
\details{
These interface to sim1fact and sim1fact, with discretization based on ucuts.
}
\examples{ 
ucuts=matrix(c(.3,.6,.4,.7,.5,.8),2,3)
param=c(5.5,6.5,4)
set.seed(123)
ydat=sim1irfact(ucuts,n=1000,param,qcond1=qcondfrk,copname1="frank")
print(cor(ydat))
for(j in 1:length(param)) print(table(ydat[,j]))
ydat2=sim2irfact(ucuts,n=1000,param,c(2,2,2),qcond1=qcondgum,qcond2=qcondgum,ivect=TRUE)
print(cor(ydat2))
for(j in 1:length(param)) print(table(ydat2[,j]))
}
\keyword{factor model}
\keyword{simulation}
