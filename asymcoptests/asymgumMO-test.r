# test of functions for asymmetric Gumbel with Marshall-Olkin at boundary

library(CopulaModel)
#source("../R/asymgumMO.R")

cat("\nchecking derivatives of A\n")
heps=1.e-5
param=c(2,.1,.2)
x=1
y=2
out=AasymgumMO(x,y,param)
outx=AasymgumMO(x+heps,y,param)
outy=AasymgumMO(x,y+heps,param)
outxy=AasymgumMO(x+heps,y+heps,param)
cat("analytic and numerical:\n")
print(out)
print((outx$Afn-out$Afn)/heps)
print((outy$Afn-out$Afn)/heps)
print((outxy$Afn-outy$Afn-outx$Afn+out$Afn)/heps^2)

cat("\nchecking pcond of asymmetric Gumbel\n")
u=.1
v=.3
#v=.8
cdf=pasymgumMO(u,v,param)
ccdf21=pcondasymgumMO21(v,u,param)
ccdf12=pcondasymgumMO12(u,v,param)
cdfu=pasymgumMO(u+heps,v,param)
cdfv=pasymgumMO(u,v+heps,param)
cat("two different conditional cdfs\n")
cat((cdfu-cdf)/heps,ccdf21,"\n")
cat((cdfv-cdf)/heps,ccdf12,"\n")

# check with chkcopderiv()
chkcopderiv(u,v,param,bcdf=pasymgumMO,pcond=pcondasymgumMO21,bpdf=dasymgumMO,
  str="asymgumMO",eps=heps)

cat("\ntau and rhoS\n")
cpar=c(2,.2,.4)
tau=asymgumMO.cpar2tau(cpar)
rho=asymgumMO.cpar2rhoS(cpar)
cat(tau,rho,"\n")

cat("\ntau and rhoS, boundary case of Gumbel\n")
cpar=c(2,1,1) 
tau=asymgumMO.cpar2tau(cpar)
rho=asymgumMO.cpar2rhoS(cpar)
cat(tau,rho,"\n")
tau=gum.cpar2tau(cpar[1])
rho=gum.cpar2rhoS(cpar[1])
cat(tau,rho,"\n")

cat("\ntau and rhoS, boundary case of independence\n")
cpar=c(2,.0,.0) # fails
cpar=c(2,.001,.001)
tau=asymgumMO.cpar2tau(cpar)
rho=asymgumMO.cpar2rhoS(cpar)
cat(tau,rho,"\n") # should be close to 0 for independence copula

