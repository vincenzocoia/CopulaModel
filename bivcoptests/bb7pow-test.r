# BB7 with power parameter

library(CopulaModel)

cat("check cop derivs for BB7\n")
param=c(1.2,1.4)
u=.3
v=seq(.4,.9,.1)
chkcopderiv(u,v,param,bcdf=pbb7,pcond=pcondbb7,bpdf=dbb7,str="bb7",eps=1.e-5)
cat("check cop derivs for BB7r\n")
param=c(1.2,1.4)
u=.6
v=seq(.4,.9,.1)
chkcopderiv(u,v,param,bcdf=pbb7r,pcond=pcondbb7r,bpdf=dbb7r,str="bb7r",eps=1.e-5)

param=c(1.3,1.4,3)
u=seq(.1,.6,.1)
v=seq(.4,.9,.1)
cat("\n")
print(pbb7pow(u,v,param))
print(pcondbb7pow(v,u,param))
print(dbb7pow(u,v,param))

cat("\ncheck cop derivs for BB7pow\n")
param=c(1.2,1.4,2.4)
#param=c(1.2,1.3,.5)
#param=c(1.3,1.4,3)
u=.3
#u=.8
v=seq(.4,.9,.1)
chkcopderiv(u,v,param,bcdf=pbb7pow,pcond=pcondbb7pow,bpdf=dbb7pow,str="bb7pow",eps=1.e-5)

