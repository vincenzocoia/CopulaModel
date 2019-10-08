# simulation from 2-factor model given a parametric bivariate copula family
# MLE for 2-factor copula log-likelihood

library(CopulaModel)

d=7
be1=c(.7,.6,.7,.6,.7,.6,.7)
be2=c(.4,.4,.4,.4,.3,.3,.3)

cpar1.frk=frk.b2cpar(be1)
cpar2.frk=frk.b2cpar(be2)
cpar1.gum=gum.b2cpar(be1)
cpar2.gum=gum.b2cpar(be2)
lmbb1=matrix(c(.3,.4,.5,.3,.5,.4,.5, .6,.6,.6,.7,.7,.6,.6),7,2)
cpar.bb1=lmbb1
for(i in 1:nrow(lmbb1))
{ cpar.bb1[i,]=bb1.lm2cpar(lmbb1[i,]) }

d=length(cpar1.gum)
n=300
set.seed(123)
frkdat=sim2fact(n,cpar1.frk,cpar2.frk,qcondfrk,qcondfrk,"frk","frk")
set.seed(123)
gumdat=sim2fact(n,cpar1.gum,cpar2.gum,qcondgum,qcondgum,"gum","gum")
set.seed(123)
bb1frkdat=sim2fact(n,cpar.bb1,cpar2.gum,qcondbb1,qcondfrk,"bb1","frk")

#cat("\nfrkdat correlations\n")
#print(cor(frkdat))
#cat("\ngumdat correlations\n")
#print(cor(gumdat))
#cat("\nbb1frkdat correlations\n")
#print(cor(bb1frkdat))

nq=31
nq=25
nq=21
cat("\nnq=",nq," for Gauss-Legendre quadrature\n")
gl = gausslegendre(nq)

# MLE
# standalone R interface and f90
cat("\nFrank 2-factor: f90, nlm and then pdhessmin\n")
out.frk=ml2fact(nq=nq,c(cpar1.frk,cpar2.frk),frkdat,copname="frank",
  LB=-30,UB=30,prlevel=1,mxiter=100)
dstrfrk=list(copname="frank",data=frkdat,quad=gl,repar=0)
out=pdhessminb(c(cpar1.frk,cpar2.frk),f90cop2nllk,ifixed=rep(F,2*d),
  dstruct=dstrfrk, LB=rep(-20,2*d),UB=rep(20,2*d),iprint=T,eps=1.e-5)
cat(cpar1.frk,cpar2.frk,"\n")
cat("\n============================================================\n")

cat("\nGumbel 2-factor: f90, nlm and then pdhessmin\n")
# some differences if nq=21 and 25 but OK for nq=31
out.gum=ml2fact(nq=nq,c(cpar1.gum,cpar2.gum),gumdat,copname="gumbel",
  LB=1,UB=20,prlevel=1,mxiter=100)
dstrgum=list(copname="gumbel",data=gumdat,quad=gl,repar=0)
out=pdhessminb(c(cpar1.gum,cpar2.gum),f90cop2nllk,ifixed=rep(F,2*d),
  dstruct=dstrgum, LB=rep(1,2*d),UB=rep(10,2*d),iprint=T,eps=1.e-5)
cat("\nstart from QN soln\n")
out=pdhessminb(out.gum$estimate,f90cop2nllk,ifixed=rep(F,2*d),
  dstruct=dstrgum, LB=rep(1,2*d),UB=rep(10,2*d),iprint=T,eps=1.e-5)
cat(cpar1.gum,cpar2.gum,"\n")
cat("\n============================================================\n")

cat("\nBB1/Frank 2-factor: f90, nlm and then pdhessmin\n")
# nq=21 not adequate,  nq=31 is OK
cpar.bb1=c(t(cpar.bb1))
out.bb1=ml2fact(nq=nq,c(cpar.bb1,cpar2.gum),bb1frkdat,copname="bb1frank",
  LB=c(rep(c(0,1),d),rep(-10,d)),UB=20,prlevel=1,mxiter=100)
dstrbb1=list(copname="bb1frank",data=bb1frkdat,quad=gl,repar=0)
out=pdhessminb(c(cpar.bb1,cpar2.frk),f90cop2nllk,ifixed=rep(F,3*d),
  dstruct=dstrbb1, LB=c(rep(c(0,1),d),rep(-10,d)),
  UB=c(rep(10,2*d),rep(20,d)),iprint=T,eps=1.e-5)
cat(cpar.bb1,cpar2.frk,"\n")
cat("\n============================================================\n")


