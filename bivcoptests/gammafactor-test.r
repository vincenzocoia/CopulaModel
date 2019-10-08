# checks of bivariate gamma factor copula 

library(CopulaModel)
#source("../R/gammaconvfactor.R")

u=.3
v=seq(.4,.9,.1)
param=c(2,1.2,1.4)
cat("\nChecking pcop, pcond, dcop\n")
chkcopderiv(u,v,param,bcdf=pbgamfcop,pcond=pcondbgamfcop21,bpdf=dbgamfcop,
  str="gammafactor",eps=1.e-5)

set.seed(123)
u=runif(1)
v=runif(5)
chkcopderiv(u,v,param,bcdf=pbgamfcop,pcond=pcondbgamfcop21,bpdf=dbgamfcop,
  str="gammafactor",eps=1.e-5)

