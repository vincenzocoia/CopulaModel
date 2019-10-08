# test of functions for bivariate Marshall-Olkin 

library(CopulaModel)
#source("../R/asymgumMO.R")

cpar=c(.1,.3)
u=.4
v=.8
vv=seq(.1,.9,.2)
uu=rep(.4,5)
out=pbMO(u,v,cpar)
print(out)
# 0.3421551
outv=pbMO(u,vv,cpar)
print(outv)
# 0.04383833 0.13151499 0.21919165 0.30686830 0.37156068
out=pcondbMO21(v,u,cpar)
print(out)
# 0.8553877
outv=pcondbMO21(vv,uu,cpar)
print(outv)
# 0.8877262 0.8877262 0.8877262 0.8877262 0.9289017
outv=pcondbMO21(vv,rep(.9,5),cpar)
print(outv)
# 0.09095326 0.27285977 0.45476629 0.63667280 0.81857932

# check if these functions are OK when one parameter is 0
cpar=c(0,.3)
outv=pbMO(u,vv,cpar)
print(outv)
# 0.04 0.12 0.20 0.28 0.36
outv=pcondbMO21(vv,uu,cpar)
print(outv)
# 0.1 0.3 0.5 0.7 0.9


