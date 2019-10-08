# checks for gppmfcdf, gp1pmfcdf, gp2pmfcdf

library(CopulaModel)
theta=5; vrh=.4
ub=10
out0=gpoispmfcdf(ub,theta,vrh)
cat("checking gpoispmfcdf\n")
pmf=dgpois(0:ub,c(theta,vrh))
cdf=rep(0,ub+1)
for(i in 0:ub) cdf[i+1]=pgpois(i,c(theta,vrh))
print(out0)
print(out0[,2]-pmf)
print(out0[,3]-cdf)

cat("\nchecking gp1pmfcdf\n")
b1=0
b0=2
x=1
vrh=0.3
xi=1/(1-vrh)^2-1
param=c(b0,b1,xi)
out1=gp1pmfcdf(ub,param,x)
mu=exp(b0+b1*x)
theta=mu*(1-vrh)
pmf=dgpois(0:ub,c(theta,vrh))
cdf=rep(0,ub+1)
for(i in 0:ub) cdf[i+1]=pgpois(i,c(theta,vrh))
print(out1)
print(out1[,2]-pmf)
print(out1[,3]-cdf)

cat("\nchecking gp2pmfcdf\n")
b1=0
b0=2
x=1
th=3.5
param=c(b0,b1,th)
out2=gp2pmfcdf(ub,param,x)
mu=exp(b0+b1*x)
vrh=1-th/mu
pmf=dgpois(0:ub,c(th,vrh))
cdf=rep(0,ub+1)
for(i in 0:ub) cdf[i+1]=pgpois(i,c(th,vrh))
print(out2)
print(out2[,2]-pmf)
print(out2[,3]-cdf)


