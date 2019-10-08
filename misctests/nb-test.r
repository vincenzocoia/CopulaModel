# checks for nbpmfcdf, nb1pmfcdf, nb2pmfcdf

library(CopulaModel)
theta=5; p=.4
ub=10
out0=nbpmfcdf(ub,theta,p)
cat("checking nbpmfcdf\n")
pmf=dnbinom(0:ub,size=theta,prob=p)
cdf=pnbinom(0:ub,size=theta,prob=p)
print(out0)
print(out0[,2]-pmf)
print(out0[,3]-cdf)

cat("\nchecking nb1pmfcdf\n")
b1=0
b0=2
xi=0.3
x=1
param=c(b0,b1,xi)
out1=nb1pmfcdf(ub,param,x)
mu=exp(b0+b1*x)
theta=mu/xi
pmf=dnbinom(0:ub,size=theta,mu=mu)
cdf=pnbinom(0:ub,size=theta,mu=mu)
print(out1)
print(out1[,2]-pmf)
print(out1[,3]-cdf)

cat("\nchecking nb2pmfcdf\n")
b1=0
b0=2
th=3.5
x=1
param=c(b0,b1,th)
out2=nb2pmfcdf(ub,param,x)
mu=exp(b0+b1*x)
pmf=dnbinom(0:ub,size=th,mu=mu)
cdf=pnbinom(0:ub,size=th,mu=mu)
print(out2)
print(out2[,2]-pmf)
print(out2[,3]-cdf)


