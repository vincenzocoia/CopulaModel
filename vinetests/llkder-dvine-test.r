# compare D-vine and R-vine codes for log-likelihood with gradient

library(CopulaModel)
#source("../R/rvinenllkderiv.R")

d=5
# choose of one following three
#qcondnames=rep("qcondfrk",4); pcondnames=rep("pcondfrk",4)
#parvec=c(3.6,3.6,3.6,3.6, 1.5,1.5,1.5, 1.4,1.4, 1.3)
#logdcopdernames=rep("logdfrk.deriv",4)
#pconddernames=rep("pcondfrk.deriv",4)
#lb=-10; ub=30

qcondnames=rep("qcondgum",4); pcondnames=rep("pcondgum",4);
parvec=c(3.6,3.6,3.6,3.6, 1.5,1.5,1.5, 1.4,1.4, 1.3)
logdcopdernames=rep("logdgum.deriv",4)
pconddernames=rep("pcondgum.deriv",4)
lb=1; ub=15

#qcondnames=rep("qcondt",4); pcondnames=rep("pcondt",4); dfdefault=5
#parvec=c(.6,.6,.6,.6, .5,.5,.5, .4,.4, .3)
#logdcopdernames=rep("logdtcop.deriv",4)
#pconddernames=rep("pcondt.deriv",4)
#lb=-1; ub=1

np=matrix(0,d,d)
np[1,2:d]=1
np[2,3:d]=1
np[3,4:d]=1
np[4,5]=1

D=Dvinearray(d)
set.seed(123)
nsim=20
udat=rvinesimvec(nsim,D,parvec,np,qcondnames,pcondnames,iprint=F)

cat("5-dimensional D-vine\n")
outd=dvinenllkder1.trunc(parvec,udat,logdcopdernames,pconddernames,LB=lb,UB=ub)
print(outd)
outr=rvinenllkder1.trunc(parvec,udat,D,logdcopdernames,pconddernames,LB=lb,UB=ub)
print(outr)
print(max(abs(outd-outr)))
print(max(abs(attr(outd,"gradient")-attr(outr,"gradient"))))

cat("\n3-truncated\n")
outd=dvinenllkder1.trunc(parvec[1:9],udat,logdcopdernames[1:3],
  pconddernames[1:3],LB=lb,UB=ub)
print(outd)
outr=rvinenllkder1.trunc(parvec[1:9],udat,D,logdcopdernames[1:3],
  pconddernames[1:3],LB=lb,UB=ub)
print(outr)
print(max(abs(outd-outr)))
print(max(abs(attr(outd,"gradient")-attr(outr,"gradient"))))

cat("\n2-truncated\n")
outd=dvinenllkder1.trunc(parvec[1:7],udat,logdcopdernames[1:2],
  pconddernames[1:2],LB=lb,UB=ub)
print(outd)
outr=rvinenllkder1.trunc(parvec[1:7],udat,D,logdcopdernames[1:2],
  pconddernames[1:2],LB=lb,UB=ub)
print(outr)
print(max(abs(outd-outr)))
print(max(abs(attr(outd,"gradient")-attr(outr,"gradient"))))

cat("\n1-truncated\n")
outd=dvinenllkder1.trunc(parvec[1:4],udat,logdcopdernames[1],
  pconddernames[1],LB=lb,UB=ub)
print(outd)
outr=rvinenllkder1.trunc(parvec[1:4],udat,D,logdcopdernames[1],
  pconddernames[1],LB=lb,UB=ub)
print(outr)
print(max(abs(outd-outr)))
print(max(abs(attr(outd,"gradient")-attr(outr,"gradient"))))

