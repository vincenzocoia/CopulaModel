# checks of bivariate GB2 copula 

library(CopulaModel)
#source("../R/genbeta2.R")

u=.3
v=seq(.4,.9,.1)
param=c(1.5,1.8,2.3)
cat("\nChecking pcop, pcond, dcop\n")
chkcopderiv(u,v,param,bcdf=pbgb2cop,pcond=pcondbgb2cop,bpdf=dbgb2cop,str="gb2",eps=1.e-5)

set.seed(123)
u=runif(1)
v=runif(5)
chkcopderiv(u,v,param,bcdf=pbgb2cop,pcond=pcondbgb2cop,bpdf=dbgb2cop,str="gb2",eps=1.e-5)

