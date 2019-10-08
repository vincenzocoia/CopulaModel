# test of functions for asymmetric Gumbel (or bilogistic)

library(CopulaModel)
#source("../R/asymgum.R")
ze=.6; eta=.6
ze=.6; eta=.7
ze=.7; eta=.6
ze=.3; eta=.2
ww=seq(.05,.95,.05)
out=Basymgum(ww,c(ze,eta))
cat("B function for ", ze, eta,"\n")
print(cbind(ww,out$Bfn,out$Bder,out$Bder2))

# check of first and second order derivatives of B
cat("\nchecking derivatives of B\n")
w=.1
heps=1.e-5
weps=w+heps
out1=Basymgum(w,c(ze,eta))
out2=Basymgum(w+heps,c(ze,eta))
out3=Basymgum(w+2*heps,c(ze,eta))
der=(out2$Bfn-out1$Bfn)/heps
der2=(out3$Bfn-2*out2$Bfn+out1$Bfn)/heps^2
cat("analytic and numerical:\n")
print(out1)
cat(der,der2,"\n")

# check pcond
cat("\nchecking pcond of asymmetric Gumbel\n")
u=.1
v=.3
#v=.8
cdf=pasymgum(u,v,c(ze,eta))
ccdf21=pcondasymgum21(v,u,c(ze,eta))
ccdf12=pcondasymgum12(u,v,c(ze,eta))
cdfu=pasymgum(u+heps,v,c(ze,eta))
cdfv=pasymgum(u,v+heps,c(ze,eta))
cat("two different conditional cdfs\n")
cat((cdfu-cdf)/heps,ccdf21,"\n")
cat((cdfv-cdf)/heps,ccdf12,"\n")
param=c(ze,eta)
chkcopderiv(u,v,param,bcdf=pasymgum,pcond=pcondasymgum21,bpdf=dasymgum,
  str="asymgum",eps=heps)

# compare Gumbel in special case
cat("\nsymmetric case Gumbel\n")
param=c(.5,.5)
cdf1=pasymgum(u,v,param)
cdf2=pgum(u,v,1/param[1])  # OK
cat(cdf1,cdf2,"\n")

v=seq(.1,.9,.2)
cdf1=pasymgum(u,v,param)
cdf2=pgum(u,v,1/param[1]) 
print(cbind(cdf1,cdf2))

cat("\ntau and rhoS, symmetric case\n")
param=c(.5,.5)
tau=asymgum.cpar2tau(param)
rho=asymgum.cpar2rhoS(param)
cat(tau,rho,"\n")
cat("\ntau and rhoS, compare Gumbel\n")
tau=gum.cpar2tau(1/param[1])
rho=gum.cpar2rhoS(1/param[1])
cat(tau,rho,"\n")

cat("\ntau and rhoS, asymmetric case\n")
#param=c(.6,.7)
param=c(.2,.4)
tau=asymgum.cpar2tau(param)
rho=asymgum.cpar2rhoS(param)
cat(tau,rho,"\n")
