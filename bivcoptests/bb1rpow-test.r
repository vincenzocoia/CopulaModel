# BB1 reflected with power parameter

library(CopulaModel)

cat("check cop derivs for BB1r\n")
param=c(1.2,1.4)
u=.3
v=seq(.4,.9,.1)
chkcopderiv(u,v,param,bcdf=pbb1r,pcond=pcondbb1r,bpdf=dbb1r,str="bb1r",eps=1.e-5)

cat("check cop derivs for BB1rpow\n")
param=c(1.2,1.4,2.4)
param=c(1.2,1.3,.5)
param=c(1.3,1.4,3)
u=.3
#u=.8
v=seq(.4,.9,.1)
chkcopderiv(u,v,param,bcdf=pbb1rpow,pcond=pcondbb1rpow,bpdf=dbb1rpow,str="bb1rpow",eps=1.e-5)

