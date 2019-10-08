# 1-factor model for item response, 
# Science data in Nikoloulopoulos and Joe (2014)

# QN using nlm, and pdhessmin
library(CopulaModel)
library(abind)
data(ltmconv)
d=ncol(sci)
nq=21

# using nlm and quasi-Newton
cat("Science data with 1-factor Gumbel, nlm with R function\n")
out=ml1irfact(nq,start=rep(2,d),sci,pcond=pcondgum,LB=1,UB=20,ihess=T,prlevel=1,mxiter=50)
theta=out$est 
tau=gum.cpar2tau(theta)
cat("\ntaus and SEs:\n")
print(tau)
# SE conversion SE(tau)=SE(theta)/theta^2
print(out$SE/theta^2)  

cat("\nScience data with 1-factor Gumbel, nlm with f90 function\n")
out=f90ml1irfact(nq,start=rep(2,d),sci,copname="gumbel",LB=1,UB=20,ihess=T,prlevel=1,mxiter=50)
# differences in 5th/6th decimal place

cat("\n============================================================\n")

cat("\nScience data with 1-factor reflected Gumbel, nlm with R function\n")
asci=3-sci
outr=ml1irfact(nq,start=rep(2,d),asci,pcond=pcondgum,LB=1,UB=20,prlevel=1,mxiter=50)
thetar=outr$est 
taur=gum.cpar2tau(thetar)
cat("\ntaus and SEs:\n")
print(taur)
print(outr$SE/thetar^2)  # nothing if not posdef

cat("\n============================================================\n")

cat("Science data with 1-factor Gaussian, nlm with R function\n")
out=ml1irfact(nq,start=rep(.5,d),sci,pcond=pcondbvncop,LB=-1,UB=1,ihess=T,prlevel=1,mxiter=50)
theta=out$est 
tau=bvn.cpar2tau(theta)
cat("\ntaus:\n")
print(tau)

cat("\nScience data with 1-factor Gaussian, nlm with f90 function\n")
out=f90ml1irfact(nq,start=rep(.5,d),sci,copname="gaussian",LB=-1,UB=1,ihess=T,prlevel=1,mxiter=50)
# differences in 6th decimal place

cat("\n============================================================\n")

cat("Science data with 1-factor t(2), nlm with R function\n")
dfdefault=2
out=ml1irfact(nq,start=rep(.5,d),sci,pcond=pcondt,LB=-1,UB=1,prlevel=1,mxiter=50)
theta=out$est 
tau=bvn.cpar2tau(theta)
cat("\ntaus:\n")
print(tau)

cat("\nScience data with 1-factor t(2), nlm with f90 function\n")
out=f90ml1irfact(nq,start=rep(.5,d),sci,copname="t",nu=2,LB=-1,UB=1,ihess=T,prlevel=1,mxiter=50)

# matching Table 6

cat("\n============================================================\n")


#============================================================

# modified Newton-Raphson and f90 code

ucutp=unifcuts(sci)
gl = gausslegendre(nq)
dstrgum=list(copname="gumbel",dat=sci,quad=gl,cutp=ucutp)  

cat("\n1-factor Gumbel: pdhessmin\n")
ir1nllk(param=rep(2,d),dstruct=dstrgum,iprfn=T)
ir1nllk(param=c(1.5,1.1,1.6,2.5,1.05,1.2,1.5),dstruct=dstrgum,iprfn=T)
outg1 = pdhessmin(param=rep(2,d),ir1nllk,dstruct=dstrgum,LB=rep(1,d),UB=rep(20,d),iprint=T,eps=1.e-5);
print(outg1$iposdef) # T
#outg1 = pdhessmin(param=c(1.5,1.1,1.6,2.5,1.05,1.2,1.5),ir1nllk,dstruct=dstrgum,LB=rep(1,d),UB=rep(20,d),iprint=T,eps=1.e-5);
#print(outg1$iposdef)
print(outg1$parmin);
print(sqrt(diag(outg1$invh)))

cat("\n============================================================\n")

mleqn=c(.48,-.02,.55,.80,-.02,.14,.50)
dstrgau=list(copname="gaussian",dat=sci,quad=gl,cutp=ucutp)  
cat("\n1-factor Gaussian: pdhessmin\n")
ir1nllk(param=rep(.5,d),dstruct=dstrgau,iprfn=T)
ir1nllk(param=mleqn,dstruct=dstrgau,iprfn=T)
outn1 = pdhessmin(param=rep(.5,d),ir1nllk,dstruct=dstrgau,LB=rep(-1,d),UB=rep(1,d),iprint=T,eps=1.e-5);
print(outn1$iposdef) # T
#outn1 = pdhessmin(param=mleqn,ir1nllk,dstruct=dstrgau,LB=rep(-1,d),UB=rep(1,d),iprint=T,eps=1.e-5);
#print(outn1$iposdef)
print(outn1$parmin);
print(sqrt(diag(outn1$invh)))

cat("\n============================================================\n")

mleqn=c(.48,-.02,.55,.80,-.02,.14,.50)
dstrt=list(copname="t",dat=sci,quad=gl,cutp=ucutp,nu=2)  
cat("\n1-factor t(2): pdhessmin\n")
outt1 = pdhessmin(param=mleqn,ir1nllk,dstruct=dstrt,LB=rep(-1,d),UB=rep(1,d),iprint=T,eps=1.e-5);
print(outt1$iposdef) # 
print(outt1$parmin);
print(sqrt(diag(outt1$invh)))

