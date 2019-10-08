# 2-factor model for item response, 
# Science data in Nikoloulopoulos and Joe (2014)
# QN using nlm, to be compared with pdhessmin

library(CopulaModel)
library(abind)
data(ltmconv)
d=ncol(sci)
ucutp=unifcuts(sci)
np=2*d
nq=21
gl = gausslegendre(nq)

# using nlm and quasi-Newton
cat("Science data with 2-factor Gumbel, R/nlm\n")
param=c(1.5,1.1,1.6,2.5,1.05,1.2,1.5,rep(1.1,d))
out=ml2irfact(nq,start=param,sci,pcond1=pcondgum,pcond2=pcondgum,LB=1,UB=20,prlevel=1,mxiter=50)

hess=out$hessian
print(isposdef(hess)) # F ( one parameter at boundary)
qnmle.gum=out$estimate

# compare modified Newton-Raphson
cat("\nScience data with 2-factor Gumbel, pdhessmin\n")
dstrgum=list(copname="gumbel",dat=sci,quad=gl,cutp=ucutp) 
print(ir2nllk(param,dstruct=dstrgum,iprfn=T))
print(ir2nllk(qnmle.gum,dstruct=dstrgum,iprfn=T)) # check if match last iter of nlm
outg2 = pdhessmin(param,ir2nllk,dstruct=dstrgum,LB=rep(1,np),UB=rep(20,np),iprint=T,eps=1.e-5);

cat("\n============================================================\n")

cat("\nScience data with 2-factor Gumbel/t2, R/nlm\n")
param=c(1.5,1.1,1.6,2.5,1.05,1.2,1.5,rep(.4,d))
dfdefault=2
out=ml2irfact(nq,start=param,sci,pcond1=pcondgum,pcond2=pcondt,
  LB=c(rep(1,d),rep(-1,d)),UB=c(rep(20,d),rep(1,d)),prlevel=1,mxiter=50)
cat("\nScience data with 2-factor Gumbel/t2, pdhessmin\n")
dstrgumt=list(copname="gumbelt",dat=sci,quad=gl,cutp=ucutp,nu2=2) 
#tem=ir2nllk(param,dstrgumt,iprfn=T)
outgt = pdhessmin(param,ir2nllk,dstruct=dstrgumt,
  LB=c(rep(1,d),rep(-1,d)),UB=c(rep(20,d),rep(1,d)),iprint=T,eps=1.e-5);
cat("\nScience data with 2-factor Gumbel/t2, f90/nlm\n")
outnlm=f90ml2irfact(nq,param,sci,"gumbelt",LB=c(rep(1,d),rep(-1,d)),
  UB=c(rep(20,d),rep(1,d)), ihess=F,prlevel=1,mxiter=50,nu2=2) 
cat("\n============================================================\n")

cat("\nScience data with 2-factor t2/Gumbel, R/nlm\n")
param=c(rep(.4,d),1.5,1.1,1.6,2.5,1.05,1.2,1.5)
dfdefault=2
out=ml2irfact(nq,start=param,sci,pcond1=pcondt,pcond2=pcondgum,
  LB=c(rep(-1,d),rep(1,d)),UB=c(rep(1,d),rep(20,d)),prlevel=1,mxiter=50)
cat("\nScience data with 2-factor t2/Gumbel, pdhessmin\n")
dstrtgum=list(copname="tgumbel",dat=sci,quad=gl,cutp=ucutp,nu1=2) 
#tem=ir2nllk(param,dstrtgum,iprfn=T)
outtg = pdhessmin(param,ir2nllk,dstruct=dstrtgum,
  LB=c(rep(-1,d),rep(1,d)),UB=c(rep(1,d),rep(20,d)),iprint=T,eps=1.e-5);
cat("\n============================================================\n")

cat("\nScience data with 2-factor t3/Gumbel, R/nlm\n")
param=c(rep(.4,d),1.5,1.1,1.6,2.5,1.05,1.2,1.5)
dfdefault=3
out=ml2irfact(nq,start=param,sci,pcond1=pcondt,pcond2=pcondgum,
  LB=c(rep(-1,d),rep(1,d)),UB=c(rep(1,d),rep(20,d)),prlevel=1,mxiter=50)
cat("\nScience data with 2-factor t3/Gumbel, pdhessmin\n")
dstrtgum=list(copname="tgumbel",dat=sci,quad=gl,cutp=ucutp,nu1=3) 
#tem=ir2nllk(param,dstrtgum,iprfn=T)
outtg = pdhessmin(param,ir2nllk,dstruct=dstrtgum,
  LB=c(rep(-1,d),rep(1,d)),UB=c(rep(1,d),rep(20,d)),iprint=T,eps=1.e-5);
cat("\n============================================================\n")

cat("\nScience data with 2-factor t5/t5, R/nlm\n")
param=c(rep(.5,d),rep(.3,d))
dfdefault=5
out=ml2irfact(nq,start=param,sci,pcond1=pcondt,pcond2=pcondt,
  LB=c(rep(-1,d),rep(-1,d)),UB=c(rep(1,d),rep(1,d)),prlevel=1,mxiter=50)
cat("\nScience data with 2-factor t5/t5, pdhessmin\n")
dstrtt=list(copname="t",dat=sci,quad=gl,cutp=ucutp,nu1=5,nu2=5) 
#tem=ir2nllk(param,dstrtt,iprfn=T) 
outtt = pdhessmin(param,ir2nllk,dstruct=dstrtt,
  LB=c(rep(-1,d),rep(-1,d)),UB=c(rep(1,d),rep(1,d)),iprint=T,eps=1.e-5);
cat("\nScience data with 2-factor t5/t5, f90/nlm\n")
outnlm=f90ml2irfact(nq,param,sci,"t",LB=c(rep(-1,d),rep(-1,d)),
  UB=c(rep(1,d),rep(1,d)), ihess=F,prlevel=1,mxiter=50,nu1=5,nu2=5) 
cat("\n============================================================\n")

#============================================================

# pdhessminb
ifixed=rep(F,np)
ifixed[10]=T
start2=qnmle.gum
start2[10]=1.02
cat("\npdhessminb: copula for variable ", 3, " and factor 2 = Cindep\n")
outg2 = pdhessminb(start2,ir2nllk,ifixed,dstruct=dstrgum,LB=rep(1,np),UB=rep(20,np),iprint=T,eps=1.e-5);
print(outg2$iposdef)  

