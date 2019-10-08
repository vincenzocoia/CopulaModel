# inf for bb2, bb3, bb7
library(CopulaModel)

# check infinite case for BB2
cat("\nBB2 infinite x  or y\n")
v=.044325
u=.044316
eps=1.e-8
param=c(3,.3)
cdf=pbb2(u,v,param,iprint=T)
cdf2=pbb2(u+eps,v,param,iprint=T)
ccdf=pcondbb2(v,u,param,iprint=T)
ccdf2=pcondbb2(v+eps,u,param,iprint=T)
pdf=dbb2(u,v,param,iprint=T)
cat(u,v,cdf,(cdf2-cdf)/eps,ccdf,"\n")
cat(u,v,ccdf,(ccdf2-ccdf)/eps,pdf,"\n")

# check infinite case for BB3
cat("\nBB3 infinite x  or y\n")
v=.004603
u=.004603
eps=1.e-8
eps=1.e-5
param=c(20,2) # unreliable
param=c(10,2) # unreliable
param=c(7,2)  # closer
param=c(6,2)  # better
param=c(3,1.5)  # getting better
param=c(2.1,1.2)  # almost OK
cdf=pbb3(u,v,param,iprint=T)
cdf2=pbb3(u+eps,v,param,iprint=T)
ccdf=pcondbb3(v,u,param,iprint=T)
ccdf2=pcondbb3(v+eps,u,param,iprint=T)
pdf=dbb3(u,v,param,iprint=T)
cat(u,v,cdf,(cdf2-cdf)/eps,ccdf,"\n")
cat(u,v,ccdf,(ccdf2-ccdf)/eps,pdf,"\n")

# check x near 0 case for BB3
# separate this into separate file
cat("\nBB3 x,y near 0\n")
v=0.9402022
u=0.936651
eps=1.e-8
param=c(20,2) # OK
param=c(10,2) # OK
cdf=pbb3(u,v,param,iprint=T)
cdf2=pbb3(u+eps,v,param,iprint=T)
ccdf=pcondbb3(v,u,param,iprint=T)
ccdf2=pcondbb3(v+eps,u,param,iprint=T)
pdf=dbb3(u,v,param,iprint=T)
cat(u,v,cdf,(cdf2-cdf)/eps,ccdf,"\n")
cat(u,v,ccdf,(ccdf2-ccdf)/eps,pdf,"\n")

# checks for the difficult cases in C runs
# BB7 th=2.5 de=0.1 p=0.705077 u=9.999998e-01
vv=qcondbb7(0.705077,9.999998e-01,c(2.5,.1))
print(vv) # 0.9999998
out=pcondbb7(vv,9.999998e-01,c(2.5,.1))
print(out) # NaN

vv=qcondbb7(0.705077,.98,c(2.5,.1))
print(vv) # 0.9817972
out=pcondbb7(vv,.98,c(2.5,.1))
print(out) # 0.705077
