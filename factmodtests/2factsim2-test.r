# simulation from 2-factor model given a parametric bivariate copula family
# MLE for 2-factor copula log-likelihood

library(CopulaModel)

d=7
be1=c(.7,.6,.7,.6,.7,.6,.7)
#be2=c(.4,.4,.4,.4,.3,.3,.3)
be2=c(.1,.1,.1,.1,.2,.2,.2)

cpar1.frk=frk.b2cpar(be1)
cpar2.frk=frk.b2cpar(be2)
cpar1.gum=gum.b2cpar(be1)
cpar2.gum=gum.b2cpar(be2)

d=length(cpar1.gum)
n=300
set.seed(123)
frkdat=sim2fact(n,cpar1.frk,cpar2.frk,qcondfrk,qcondfrk,"frk","frk")
set.seed(123)
gumdat=sim2fact(n,cpar1.gum,cpar2.gum,qcondgum,qcondgum,"gum","gum")

cat("\nfrkdat correlations\n")
print(cor(frkdat))
cat("\ngumdat correlations\n")
print(cor(gumdat))

gl21 = gausslegendre(21)

# MLE: standalone R interface and f90
cat("\nFrank 2-factor: f90, nlm and then pdhessmin\n")
out.frk=ml2fact(nq=21,c(cpar1.frk,cpar2.frk),frkdat,copname="frank",
  LB=-30,UB=30,prlevel=1,mxiter=100)
dstrfrk=list(copname="frank",data=frkdat,quad=gl21,repar=0)
out=pdhessminb(c(cpar1.frk,cpar2.frk),f90cop2nllk,ifixed=rep(F,2*d),
  dstruct=dstrfrk, LB=rep(-20,2*d),UB=rep(20,2*d),iprint=T,eps=1.e-5)
cat(cpar1.frk,cpar2.frk,"\n")
cat("\n============================================================\n")

cat("\nGumbel 2-factor: f90, nlm and then pdhessmin\n")
out.gum=ml2fact(nq=21,c(cpar1.gum,cpar2.gum),gumdat,copname="gumbel",
  LB=1,UB=20,prlevel=1,mxiter=100)
dstrgum=list(copname="gumbel",data=gumdat,quad=gl21,repar=0)
out=pdhessminb(c(cpar1.gum,cpar2.gum),f90cop2nllk,ifixed=rep(F,2*d),
  dstruct=dstrgum, LB=rep(1,2*d),UB=rep(10,2*d),iprint=T,eps=1.e-5)
cat("\nstart from QN soln\n")
out=pdhessminb(out.gum$estimate,f90cop2nllk,ifixed=rep(F,2*d),
  dstruct=dstrgum, LB=rep(1,2*d),UB=rep(10,2*d),iprint=T,eps=1.e-5)
cat(cpar1.gum,cpar2.gum,"\n")
cat("\n============================================================\n")


