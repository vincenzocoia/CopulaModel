\name{vinemargincdf}
\Rdversion{1.1}
\alias{ptree2cop}
\alias{ptree3cop.cvine}
\alias{ptree3cop.dvine}
\alias{pcondtree2}
\alias{pcondtree3.cvine}
\alias{pcondtree3.dvine}
\alias{dtree2cop}
\alias{dtree3cop.cvine}
\alias{dtree3cop.dvine}
\title{
Bivariate marginal copulas for trees 2 and 3 of vine copulas
}
\description{
Bivariate marginal copulas for trees 2 and 3 of vine copulas
}
\usage{
ptree2cop(ub,uc,param,pcondba,pcondca,pcopbc,nq)
dtree2cop(ub,uc,param,pcondba,pcondca,dcopbc,dcopab,dcopac,nq)
pcondtree2(uc,ub,param,pcondba,pcondca,pcondcb,dcopab,nq)
ptree3cop.cvine(uc,ud,param,pcondba,pcondca,pcondda,pcondcb,pconddb,
  pcopcd,dcopab,nq)
pcondtree3.cvine(ud,uc,param,pcondba,pcondca,pcondda,pcondcb,pconddb,
  pconddc,dcopab,dcopac,dcopbc,nq)
dtree3cop.cvine(uc,ud,param,pcondba,pcondca,pcondda,pcondcb,pconddb,
  dcopcd,dcopab,dcopac,dcopad,dcopbc,dcopbd,nq)
ptree3cop.dvine(ua,ud,param,pcondab,pcondbc,pcondcb,pconddc,pcondac,
  pconddb,pcopad,dcopbc,nq)
pcondtree3.dvine(ud,ua,param,pcondab,pcondbc,pcondcb,pconddc,pcondac,
  pconddb,pcondda,dcopab,dcopbc,dcopac,nq)
dtree3cop.dvine(ua,ud,param,pcondab,pcondbc,pcondcb,pconddc,pcondac,
  pconddb,dcopad,dcopab,dcopbc,dcopcd,dcopac,dcopbd,nq) 
}
\arguments{
\item{ub}{vector for variable b, each element in interval 0,1}
\item{param}{vector of length 3 for tree2, and length 6 for tree3}
\item{pcopbc}{function for copula cdf for variables b,c;
similar for pcopcd, pcopad}
\item{pcondba}{function for conditional cdf of variable b given variable a;
similar for pcondca, pcondcb, pcondda, pconddb, pconddc, pcondab, pcondac,
pcondbc}
\item{dcopbc}{function for copula pdf of variables b,c;
similar for dcopab, dcopac, dcopad, dcopbd, dcopcd}
\item{nq}{number of quadrature points for Gauss-Legendre quadrature}
\item{ua}{vector for variable a, each element in interval 0,1}
\item{uc}{vector for variable c, each element in interval 0,1}
\item{ud}{vector for variable d, each element in interval 0,1}
\item{pcopcd}{similar}
\item{pcopad}{similar}
\item{pcondca}{similar}
\item{pcondcb}{similar}
\item{pcondda}{similar}
\item{pconddb}{similar}
\item{pconddc}{similar}
\item{pcondab}{similar}
\item{pcondac}{similar}
\item{pcondbc}{similar}
\item{dcopab}{similar}
\item{dcopac}{similar}
\item{dcopad}{similar}
\item{dcopbd}{similar}
\item{dcopcd}{similar}
}
\value{
cdf or conditional cdf or pdf value(s) of bivariate margin of vine copula
}
\seealso{
\code{\link{rvinenllk}}
\code{\link{pcond}}
}
\examples{
ptree2gumfrk=function(u,v,param)
{ ptree2cop(u,v,param,pcondgum,pcondgum,pfrk,nq=35) }
pcondtree2gumfrk21=function(v,u,param)
{ pcondtree2(v,u,param,pcondgum,pcondgum,pcondfrk,dgum,nq=35) }
pcondtree2gumfrk12=function(v,u,param)
{ pcondtree2(v,u,param[c(2,1,3)],pcondgum,pcondgum,pcondfrk,dgum,nq=35) }
d=4
parmat=matrix(c(0,1.4,1.4,1.6, 0,0,2.6,3.4, 0,0,0,2, 0,0,0,0),4,4,byrow=TRUE)
A=Cvinearray(d)
cat("theoretical rhoS and tau for tree 2\n")
for(j in 3:d)
{ j0=A[1,j]
  j1=A[2,j]; j2=A[j,j]
  jmin=min(j0,j1); jmax=max(j0,j1)
  param3=c(parmat[1,jmax],parmat[1,j],parmat[2,j])
  rho=rhoS(param3,ptree2gumfrk)
  tau=ktau(param3,pcond12=pcondtree2gumfrk12,pcond21=pcondtree2gumfrk21)
  cat(j,rho,tau,"\n")
}
# similar for the one pair in tree 3
}
\keyword{copula}
