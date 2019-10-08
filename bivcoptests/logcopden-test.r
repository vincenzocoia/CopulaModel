# check versions of log copula densities
library(CopulaModel)

cpar=2
uvec=seq(.1,.9,.1)

cat("\nPlackett\n")
print(log(dpla(uvec,uvec,cpar)))
print(logdpla(uvec,uvec,cpar))

cat("\nFrank\n")
print(log(dfrk(uvec,uvec,cpar)))
print(logdfrk(uvec,uvec,cpar))

cat("\nMTCJ\n")
print(log(dmtcj(uvec,uvec,cpar)))
print(logdmtcj(uvec,uvec,cpar))

cat("\nMTCJ reflected\n")
print(log(dmtcj(1-uvec,1-uvec,cpar)))
print(logdmtcjr(uvec,uvec,cpar))

cat("\nJoe/B5\n")
print(log(djoe(uvec,uvec,cpar)))
print(logdjoe(uvec,uvec,cpar))

cat("\nGumbel\n")
print(log(dgum(uvec,uvec,cpar)))
print(logdgum(uvec,uvec,cpar))

cat("\nGalambos\n")
print(log(dgal(uvec,uvec,cpar)))
print(logdgal(uvec,uvec,cpar))

cat("\nHuelser-Reiss\n")
print(log(dhr(uvec,uvec,cpar)))
print(logdhr(uvec,uvec,cpar))

# tests for logdbvncop logdbvtcop
cat("\nbvn cop\n")
rh=.6
print(log(dbvncop(uvec,uvec,rh)))
print(logdbvncop(uvec,uvec,rh))

cat("\nbivariate Gaussian\n")
zvec=seq(-3,3,.5)
print(log(dbvn2(zvec,zvec,rh)))
print(logdbvn(zvec,zvec,rh))

cat("\nbvt cop\n")
nu=3.4; rh=.6
print(log(dbvtcop(uvec,uvec,c(rh,nu))))
print(logdbvtcop(uvec,uvec,c(rh,nu)))

cat("\nBB1\n")
th=2; de=2 
print(log(dbb1(uvec,uvec,c(th,de))))
print(logdbb1(uvec,uvec,c(th,de)))

