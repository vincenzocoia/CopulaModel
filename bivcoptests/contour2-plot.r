# contour densities plots with N(0,1) margins and BB1-BB3
library(CopulaModel)

par(mfrow=c(2,2),oma=c(.5,.5,.5,.5))
zvec=seq(-3,3,.2)
param=matrix(0,4,2)

param[1,]=bb1.lm2cpar(c(.4,.6))
param[2,]=bb1.lm2cpar(c(.6,.4))
param[3,]=bb1.lm2cpar(c(.3,.7))
param[4,]=bb1.lm2cpar(c(.7,.3))
for(i in 1:nrow(param))
{ contourBivCop(param[i,],zvec,dcop=dbb1)
  lm=bb1.cpar2lm(param[i,])
  title(paste("lml=",round(lm[1],2), " lmu=",round(lm[2],2)))
}
mtext("BB1 copula density with N(0,1) margins",side=1,outer=T,line=-2)

param[1,]=c(.2,.6)
param[2,]=c(.2,1.2)
param[3,]=c(.5,.6)
param[4,]=c(.5,1.3)  # maybe bimodal
for(i in 1:nrow(param))
{ contourBivCop(param[i,],zvec,dcop=dbb2)
  title(paste("th=",round(param[i,1],2), " de=",round(param[i,2],2)))
}
mtext("BB2 copula density with N(0,1) margins",side=1,outer=T,line=-2)

param[1,]=c(1.8,.6)   # maybe bimodal
param[2,]=c(1.8,1.2)
param[3,]=c(2.5,.6)
param[4,]=c(2.5,1.3)  
for(i in 1:nrow(param))
{ contourBivCop(param[i,],zvec,dcop=dbb3)
  title(paste("th=",round(param[i,1],2), " de=",round(param[i,2],2)))
}
mtext("BB3 copula density with N(0,1) margins",side=1,outer=T,line=-2)

