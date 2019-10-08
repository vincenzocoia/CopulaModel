# checks for ipsA (integrated positive stable LT for Archimedean)

library(CopulaModel)

cat("checking Kendall tau\n")
devec=c(.1,.5,.7,1,2,3)
cat("cpar=", devec,"\n")
tau=ipsA.cpar2tau(devec)
cat(tau,"\n")

cat("tau=seq(-.9,.9,.1)\n")
cpar=ipsA.tau2cpar(seq(-.9,.9,.1))
cat("cpar=", cpar,"\n")
cat("apply tau function: ", ipsA.cpar2tau(cpar),"\n")

# integral for Kendall tau formula
# s = positive value
# param = parameter >0
psiderfn=function(s,param)
{ der=exp(-s^(1/param))/gamma(1+de)
  s*der^2
}

cat("tau from 1-dimensional integral\n")
for(de in devec)
{ out=integrate(psiderfn,0,Inf,param=de)
  print(c(de,1-4*out$value))  # OK matches
}
  
cat("\nchecking pcond/qcond\n")
u=seq(.1,.6,.1)
v=seq(.4,.9,.1)
de=.6
pp=pcondipsA(v,u,de)
vv=qcondipsA(pp,u,de)
print(cbind(u,v,pp,vv))
p=v
vv=qcondipsA(p,u,de)
pp=pcondipsA(vv,u,de)
print(cbind(u,p,vv,pp))

cat("\nchecking pcop, pcond, dcop\n")
de=2
de=1.6
u=.3
v=seq(.4,.9,.1)
chkcopderiv(u,v,de,bcdf=pipsA,pcond=pcondipsA,bpdf=dipsA,str="ipsA",eps=1.e-5)


cat("\nchecking rng, nscore and semicor\n")
# nscores of random pairs
#par(mfrow=c(2,2),oma=c(.5,.5,.5,.5))
n=3000
param=2
#param=10
#param=.4
#param=.1
set.seed(123); udata=ripsA(n,cpar=param)
zz=nscore(udata)
#plot(udata); title("cpar=2"); plot(zz)
print(semicor(zz))
param=6
set.seed(123); udata=ripsA(n,cpar=param)
zz=nscore(udata)
#plot(udata); title("cpar=6"); plot(zz)
print(semicor(zz))

cat("\nKL divergence versus bivariate Gaussian\n")
tauv=seq(-.9,.9,.1)
th.bvn=bvn.b2cpar(tauv)
for(i in c(17,16,15,14,13,12,8,7,6,5,4))
{ cat("tau=",tauv[i],"\n")
  kl=KLcopvsbvn(rh=th.bvn[i],dcop2=dfrk,param2=cpar[i],copname2="ipsA",UB=7,iprint=T)
}
