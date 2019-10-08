# tests for functions for copula parameter to dependence measure and vice versa
library(CopulaModel)

dm=seq(.1,.9,.1)

# Plackett
#cpar=((1+b)/(1-b))^2
#print(d)
# 1.493827   2.250000   3.448980   5.444444   9.000000  16.000000  32.111111
# 81.000000 361.000000

cat("\nPlackett beta and rhoS\n")
cpar=pla.b2cpar(dm)
beta=4*ppla(.5,.5,cpar)-1
print(rbind(cpar,beta))
cpar=pla.rhoS2cpar(dm)
rho=pla.cpar2rhoS(cpar)
print(rbind(cpar,rho))

cat("\nFrank beta and tau\n")
cpar=frk.b2cpar(dm)
beta=4*pfrk(.5,.5,cpar)-1
nn=length(cpar); tau=rep(0,nn)
for(i in 1:nn) { tau[i]=frk.cpar2tau(cpar[i]) }
print(rbind(cpar,beta,tau))

cat("\nMTCJ beta, tau and lambda\n")
cpar=mtcj.b2cpar(dm)
beta=4*pmtcj(.5,.5,cpar)-1
print(rbind(cpar,beta))
cpar=mtcj.tau2cpar(dm)
tau=mtcj.cpar2tau(cpar)
print(rbind(cpar,tau,mtcj.cpar2lm(cpar)))

cat("\nJoe/B5 beta and tau\n")
cpar=joe.b2cpar(dm)
beta=4*pjoe(.5,.5,cpar)-1
print(rbind(cpar,beta,frk.cpar2tau(cpar)))

# Gumbel
#ln2=log(2)
#cpar=ln2/log(2-(log(1+be))/ln2)
#print(cpar)
#  1.114532 1.255384 1.434065 1.669696 1.996644 2.483585 3.290705 4.898489
# 9.709230

cat("\nGumbel bivariate, tau, rhoS, lambda\n")
cpar=gum.b2cpar(dm)
beta=4*pgum(.5,.5,cpar)-1
print(rbind(cpar,beta))
cpar=gum.tau2cpar(dm)
tau=gum.cpar2tau(cpar)
nn=length(cpar); rho=rep(0,nn)
for(i in 1:nn) { rho[i]=gum.cpar2rhoS(cpar[i]) }
print(rbind(cpar,tau,rho,gum.cpar2lm(cpar)))

# Galambos
#cpar=ln2/log(ln2/log(1+be))
#print(de)
# 0.3493499 0.5190285 0.7134752 0.9590723 1.2926845 1.7841537 2.5943013
# 4.2039930 9.0157662

cat("\nGalambos beta, tau, rhoS, lambda\n")
cpar=gal.b2cpar(dm)
beta=4*pgal(.5,.5,cpar)-1
print(rbind(cpar,beta))
nn=length(cpar); rho=rep(0,nn); tau=rep(0,nn)
for(i in 1:nn) 
{ rho[i]=gal.cpar2rhoS(cpar[i]) 
  tau[i]=gal.cpar2tau(cpar[i])
}
print(rbind(cpar,tau,rho,gal.cpar2lm(cpar)))

# H-R
#cpar=1/qnorm(1-(log(1+be))/(2*ln2))
#print(cpar)
# 0.6733317  0.8934571  1.1355352  1.4334844  1.8309972  2.4090835  3.3532216
# 5.2171685 10.7666570

cat("\nHuesler-Reiss beta, tau, rhoS and lambda\n")
cpar=hr.b2cpar(dm)
beta=4*phr(.5,.5,cpar)-1
print(rbind(cpar,beta))
nn=length(cpar); rho=rep(0,nn); tau=rep(0,nn)
for(i in 1:nn) 
{ rho[i]=hr.cpar2rhoS(cpar[i]) 
  tau[i]=hr.cpar2tau(cpar[i])
}
print(rbind(cpar,tau,rho,hr.cpar2lm(cpar)))
lm=seq(.1,.9,.1)
cpar=hr.lm2cpar(lm)
tdp=hr.cpar2lm(cpar)
print(rbind(cpar,tdp))

