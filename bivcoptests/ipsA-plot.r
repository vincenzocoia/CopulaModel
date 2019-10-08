# plots for ipsA (integrated positive stable LT for Archimedean)

library(CopulaModel)


cat("\nchecking rng, nscore and semicor\n")
# nscores of random pairs
par(mfrow=c(2,2),oma=c(.5,.5,.5,.5))
n=3000
param=2
#param=10
#param=.4
#param=.1
set.seed(123); udata=ripsA(n,cpar=param)
zz=nscore(udata)
plot(udata); title("cpar=2"); plot(zz)
print(semicor(zz))
param=6
set.seed(123); udata=ripsA(n,cpar=param)
zz=nscore(udata)
plot(udata); title("cpar=6"); plot(zz)
print(semicor(zz))

# contour density plots
zvec=seq(-3,3,.2)
devec=c(.3,.5,2,6)

par(mfrow=c(2,2),oma=c(.5,.5,.5,.5))
for(i in 1:length(devec))
{ contourBivCop(devec[i],zvec,dcop=dipsA)
  tau=ipsA.cpar2tau(devec[i])
  title(paste("de=",round(devec[i],3)," tau=",round(tau,3)))
}
mtext("ipsA copula density with N(0,1) margins",side=1,outer=T,line=-2)

