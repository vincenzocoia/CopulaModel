# check of contour density plots with N(0,1) margins
library(CopulaModel)

zvec=seq(-3,3,.2)
dm=seq(.1,.9,.1)
# sequence of tau or rhoS values

par(mfrow=c(3,3),oma=c(.5,.5,.5,.5))

for(i in 1:length(dm))
{ cpar=depmeas2cpar(dm[i],"rhoS","frank")
  contourBivCop(cpar,zvec,dcop=dfrk)
  title(paste("cpar=",round(cpar,3)," rhoS=",round(dm[i],2)))
}
mtext("Frank copula density with N(0,1) margins",side=1,outer=T,line=-2)

for(i in 1:length(dm))
{ cpar=depmeas2cpar(dm[i],"rhoS","gumbel")
  contourBivCop(cpar,zvec,dcop=dgum)
  title(paste("cpar=",round(cpar,3)," rhoS=",round(dm[i],2)))
}
mtext("Gumbel copula density with N(0,1) margins",side=1,outer=T,line=-2)

for(i in 1:length(dm))
{ cpar=depmeas2cpar(dm[i],"rhoS","joe")
  contourBivCop(cpar,zvec,dcop=djoe)
  title(paste("cpar=",round(cpar,3)," rhoS=",round(dm[i],2)))
}
mtext("Joe/B5 copula density with N(0,1) margins",side=1,outer=T,line=-2)

# plots with BB1 (separate file)
