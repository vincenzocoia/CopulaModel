# simulation from bi-factor model 
# MLE for bi-factor copula log-likelihood

library(CopulaModel)
#source("../R/structcopsim.R")

options(digits=5)
grsize=c(4,4)
d=sum(grsize)
mgrp=length(grsize)
# number of params is 2* sum(grsize)
n=600 
#n=1000
th1=frk.b2cpar(.5) # 4.875
th2=frk.b2cpar(.6) # 6.5626
th2=frk.b2cpar(.4) # 3.6021
rh1=bvn.b2cpar(.5) # 0.7071
rh2=bvn.b2cpar(.4) # 0.58779

rr=bifct(grsize,rep(rh1,d),rep(rh2,d))
print(rr$fctmat)  # .67  within block, .5 between blocks

set.seed(123)
# Frank for observed with global, and with group latent variable
udat=simbifact(n, grsize, cop=5, c(rep(th2,d),rep(th1,d)))
r=cor(udat)
print(r) # .4-.6 within block, .2-.4 between blocks

cat("\nFrank bi-factor\n")
gl=gausslegendre(25)
dstrfrk=list(data=udat,copname="frank",quad=gl,repar=0,grsize=grsize,pdf=0);
npar=2*d
param=rep(3,npar)
LB=rep(-10,npar)
UB=rep(30,npar)
ifixed=rep(F,npar);
out=pdhessminb(param,f90str2nllk, ifixed=ifixed, dstrfrk, LB, UB, mxiter=30, 
  eps=5.e-5,iprint=T)
print(out$iposdef)
cat("SEs:\n")
print(sqrt(diag(out$invh)))
cat("parameters for simulation:\n", c(rep(th2,d),rep(th1,d)),"\n")
cat("\n============================================================\n")

# add example with multivariate t and bi-factor correlation structure
