# checks for discrete vine with bivariate Gaussian pair-copulas

library(CopulaModel)
d=4
A=Dvinearray(d)
out=varray2M(A); M=out$mxarray
ucuts4=matrix(c(.4,.5,.4,.3, .7,.8,.6,.6),2,4,byrow=T)
ucuts=rbind(rep(0.00001,d),ucuts4,rep(.99999,d))
parvec=c(.5,.5,.5,.1,.1,0)
pr=rep(1/81,81)
out=f90rvineKL(parvec, ucuts, A, M, pr)
print(out)
# 0.3789605 

data(ordinalex)
xvec=c(t(ordinalex$xx))
yvec=c(t(ordinalex$yy))
out=ordprobit.univar(xvec,yvec,iprint=F)
latentdat=mord2uu(xvec,yvec,nrep=4,out$cutpts,out$beta)
uudat=latentdat$uudat
zzdat=latentdat$zzdat

C4=Cvinearray(4)
D4=Dvinearray(4)
param=c(.5,.25,.5,.125,.25,.5)
out=varray2M(C4); MC4=out$mxarray
tem1z=rvinediscbvnnllk(param,zzdat,C4)
print(tem1z)
# 847.9132 
out=varray2M(D4); MD4=out$mxarray
tem2z=rvinediscbvnnllk(param,zzdat,D4)
print(tem2z)
# 853.781
mlz1=nlm(rvinediscbvnnllk,p=param,zzdat=zzdat,A=C4,hessian=T,print.level=1)
#iteration = 8
#Parameter:
#[1] 0.36887825 0.16064335 0.02171784 0.44478355 0.34654634 0.35107517
#Function Value
#[1] 819.1659
#Gradient:
#[1]  2.169145e-04  5.229595e-05 -4.222329e-04  3.387868e-05 -1.826947e-04
#[6]  3.213927e-04

fparam=c(-0.5285182,0.5041516,0.3642926,
   0.37107109,0.46115563,0.43073832,0.03254656,0.20708160,-0.10396673)
fnllk=rvinediscbvnfullnllk(fparam,D4,xvec,yvec,nrep=4,ncateg=3)
print(fnllk)
# 818.8315
fmlz=nlm(rvinediscbvnfullnllk,p=fparam,xmat=xvec,yvec=yvec,A=D4,
  nrep=4,ncateg=3, hessian=T,print.level=1)
#iteration = 9
#Parameter:
#[1] -0.51766478  0.51393623  0.39621155  0.37171336  0.46012778  0.43498099
#[7]  0.03577443  0.20344139 -0.10423191
#Function Value
#[1] 818.6813
#Gradient:
#[1]  1.733724e-04 -3.198011e-04  1.608669e-04 -6.298251e-05  2.953584e-04
#[6]  4.683898e-05 -1.675744e-04 -2.332854e-04  1.938361e-04

hess=fmlz$hess
acov=solve(hess)
cat("SEs\n")
print(sqrt(diag(acov)))
#[1] 0.05762959 0.05737333 0.06379935 0.08012361 0.07086124 0.07365039 0.09240526
#[8] 0.08847101 0.09207092

# test example with second random covariate
set.seed(123)
xvec2=runif(800,-1,1)
xmat=cbind(xvec,xvec2)
fparam2=c(-0.5285182,0.5041516,0.3642926,0.1,
   0.37107109,0.46115563,0.43073832,0.03254656,0.20708160,-0.10396673)
fnllk2=rvinediscbvnfullnllk(fparam2,D4,xmat,yvec,nrep=4,ncateg=3)
print(fnllk2)
# 818.649
fmlz2=nlm(rvinediscbvnfullnllk,p=fparam2,xmat=xmat,yvec=yvec,A=D4,
  nrep=4,ncateg=3, hessian=T,print.level=1)
#iteration = 9
#Parameter:
# [1] -0.51758844  0.51461501  0.39734373  0.05966108  0.37206112  0.45438427
# [7]  0.43613973  0.03335408  0.21310673 -0.11407821
#Function Value
#[1] 818.2311
#Gradient:
# [1] -5.309175e-05 -1.242597e-04  5.036327e-05 -1.867875e-04  9.060841e-05
# [6]  1.548415e-04  2.330580e-05  1.496119e-04 -8.014922e-05  6.730261e-05

